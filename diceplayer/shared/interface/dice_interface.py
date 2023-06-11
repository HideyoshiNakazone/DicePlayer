from __future__ import annotations

from diceplayer import logger
from diceplayer.shared.config.player_config import PlayerConfig
from diceplayer.shared.environment.system import System
from diceplayer.shared.interface import Interface

from setproctitle import setproctitle

import os
import random
import shutil
import subprocess
import sys
import time
from multiprocessing import Process, connection
from pathlib import Path
from typing import Final, TextIO

DICE_END_FLAG: Final[str] = "End of simulation"
DICE_FLAG_LINE: Final[int] = -2
UMAANG3_TO_GCM3: Final[float] = 1.6605

MAX_SEED: Final[int] = 4294967295


class DiceInterface(Interface):
    title = "Diceplayer run"

    def configure(self, step: PlayerConfig, system: System):
        self.step = step
        self.system = system

    def start(self, cycle: int):
        procs = []
        sentinels = []

        for proc in range(1, self.step.nprocs + 1):
            p = Process(target=self._simulation_process, args=(cycle, proc))
            p.start()

            procs.append(p)
            sentinels.append(p.sentinel)

        while procs:
            finished = connection.wait(sentinels)
            for proc_sentinel in finished:
                i = sentinels.index(proc_sentinel)
                status = procs[i].exitcode
                procs.pop(i)
                sentinels.pop(i)
                if status != 0:
                    for p in procs:
                        p.terminate()
                    sys.exit(status)

        logger.info("\n")

    def reset(self):
        del self.step
        del self.system

    def _simulation_process(self, cycle: int, proc: int):
        setproctitle(f"diceplayer-step{cycle:0d}-p{proc:0d}")

        try:
            self._make_proc_dir(cycle, proc)
            self._make_dice_inputs(cycle, proc)
            self._run_dice(cycle, proc)
        except Exception as err:
            sys.exit(err)

    def _make_proc_dir(self, cycle, proc):
        simulation_dir = Path(self.step.simulation_dir)
        if not simulation_dir.exists():
            simulation_dir.mkdir(parents=True)

        proc_dir = Path(simulation_dir, f"step{cycle:02d}", f"p{proc:02d}")
        proc_dir.mkdir(parents=True, exist_ok=True)

    def _make_dice_inputs(self, cycle, proc):
        proc_dir = Path(self.step.simulation_dir, f"step{cycle:02d}", f"p{proc:02d}")

        self._make_potentials(proc_dir)

        random.seed(self._make_dice_seed())

        # This is logic is used to make the initial configuration file
        # for the next cycle using the last.xyz file from the previous cycle.
        if self.step.dice.randominit == "first" and cycle > 1:
            last_xyz = Path(
                self.step.simulation_dir,
                f"step{(cycle - 1):02d}",
                f"p{proc:02d}",
                "last.xyz",
            )
            if not last_xyz.exists():
                raise FileNotFoundError(f"File {last_xyz} not found.")

            with open(last_xyz, "r") as last_xyz_file:
                self._make_init_file(proc_dir, last_xyz_file)
                last_xyz_file.seek(0)
                self.step.dice.dens = self._new_density(last_xyz_file)

        else:
            self._make_nvt_ter(cycle, proc_dir)

        if len(self.step.dice.nstep) == 2:
            self._make_nvt_eq(cycle, proc_dir)

        elif len(self.step.dice.nstep) == 3:
            self._make_npt_ter(cycle, proc_dir)
            self._make_npt_eq(proc_dir)

    def _run_dice(self, cycle: int, proc: int):
        working_dir = os.getcwd()

        proc_dir = Path(self.step.simulation_dir, f"step{cycle:02d}", f"p{proc:02d}")

        logger.info(
            f"Simulation process {str(proc_dir)} initiated with pid {os.getpid()}"
        )

        os.chdir(proc_dir)

        if not (self.step.dice.randominit == "first" and cycle > 1):
            self.run_dice_file(cycle, proc, "NVT.ter")

        if len(self.step.dice.nstep) == 2:
            self.run_dice_file(cycle, proc, "NVT.eq")

        elif len(self.step.dice.nstep) == 3:
            self.run_dice_file(cycle, proc, "NPT.ter")
            self.run_dice_file(cycle, proc, "NPT.eq")

        os.chdir(working_dir)

        xyz_file = Path(proc_dir, "phb.xyz")
        last_xyz_file = Path(proc_dir, "last.xyz")

        if xyz_file.exists():
            shutil.copy(xyz_file, last_xyz_file)
        else:
            raise FileNotFoundError(f"File {xyz_file} not found.")

    @staticmethod
    def _make_dice_seed() -> int:
        num = time.time()
        num = (num - int(num)) * 1e6

        num = int((num - int(num)) * 1e6)

        return (os.getpid() * num) % (MAX_SEED + 1)

    def _make_init_file(self, proc_dir: Path, last_xyz_file: TextIO):
        xyz_lines = last_xyz_file.readlines()

        SECONDARY_MOLECULE_LENGTH = 0
        for i in range(1, len(self.step.dice.nmol)):
            SECONDARY_MOLECULE_LENGTH += self.step.dice.nmol[i] * len(
                self.system.molecule[i].atom
            )

        xyz_lines = xyz_lines[-SECONDARY_MOLECULE_LENGTH:]

        input_file = Path(proc_dir, self.step.dice.outname + ".xy")
        with open(input_file, "w") as f:
            for atom in self.system.molecule[0].atom:
                f.write(f"{atom.rx:>10.6f}  {atom.ry:>10.6f}  {atom.rz:>10.6f}\n")

            for line in xyz_lines:
                atom = line.split()
                rx = float(atom[1])
                ry = float(atom[2])
                rz = float(atom[3])
                f.write(f"{rx:>10.6f}  {ry:>10.6f}  {rz:>10.6f}\n")

            f.write("$end")

    def _new_density(self, last_xyz_file: TextIO):
        last_xyz_lines = last_xyz_file.readlines()

        box = last_xyz_lines[1].split()
        volume = float(box[-3]) * float(box[-2]) * float(box[-1])

        total_mass = 0
        for i in range(len(self.system.molecule)):
            total_mass += self.system.molecule[i].total_mass * self.step.dice.nmol[i]

        density = (total_mass / volume) * UMAANG3_TO_GCM3

        return density

    def _make_nvt_ter(self, cycle, proc_dir):
        file = Path(proc_dir, "NVT.ter")
        with open(file, "w") as f:
            f.write(f"title = {self.title} - NVT Thermalization\n")
            f.write(f"ncores = {self.step.ncores}\n")
            f.write(f"ljname = {self.step.dice.ljname}\n")
            f.write(f"outname = {self.step.dice.outname}\n")

            mol_string = " ".join(str(x) for x in self.step.dice.nmol)
            f.write(f"nmol = {mol_string}\n")

            f.write(f"dens = {self.step.dice.dens}\n")
            f.write(f"temp = {self.step.dice.temp}\n")

            if self.step.dice.randominit == "first" and cycle > 1:
                f.write(f"init = yesreadxyz\n")
                f.write(f"nstep = {self.step.altsteps}\n")
            else:
                f.write(f"init = yes\n")
                f.write(f"nstep = {self.step.dice.nstep[0]}\n")

            f.write("vstep = 0\n")
            f.write("mstop = 1\n")
            f.write("accum = no\n")
            f.write("iprint = 1\n")
            f.write("isave = 0\n")
            f.write("irdf = 0\n")

            seed = int(1e6 * random.random())
            f.write(f"seed = {seed}\n")
            f.write(f"upbuf = {self.step.dice.upbuf}")

    def _make_nvt_eq(self, cycle, proc_dir):
        file = Path(proc_dir, "NVT.eq")
        with open(file, "w") as f:
            f.write(f"title = {self.title} - NVT Production\n")
            f.write(f"ncores = {self.step.ncores}\n")
            f.write(f"ljname = {self.step.dice.ljname}\n")
            f.write(f"outname = {self.step.dice.outname}\n")

            mol_string = " ".join(str(x) for x in self.step.dice.nmol)
            f.write(f"nmol = {mol_string}\n")

            f.write(f"dens = {self.step.dice.dens}\n")
            f.write(f"temp = {self.step.dice.temp}\n")

            if self.step.dice.randominit == "first" and cycle > 1:
                f.write("init = yesreadxyz\n")
            else:
                f.write("init = no\n")

            f.write(f"nstep = {self.step.dice.nstep[1]}\n")

            f.write("vstep = 0\n")
            f.write("mstop = 1\n")
            f.write("accum = no\n")
            f.write("iprint = 1\n")

            f.write(f"isave = {self.step.dice.isave}\n")
            f.write(f"irdf = {10 * self.step.nprocs}\n")

            seed = int(1e6 * random.random())
            f.write("seed = {}\n".format(seed))

    def _make_npt_ter(self, cycle, proc_dir):
        file = Path(proc_dir, "NPT.ter")
        with open(file, "w") as f:
            f.write(f"title = {self.title} - NPT Thermalization\n")
            f.write(f"ncores = {self.step.ncores}\n")
            f.write(f"ljname = {self.step.dice.ljname}\n")
            f.write(f"outname = {self.step.dice.outname}\n")

            mol_string = " ".join(str(x) for x in self.step.dice.nmol)
            f.write(f"nmol = {mol_string}\n")

            f.write(f"press = {self.step.dice.press}\n")
            f.write(f"temp = {self.step.dice.temp}\n")

            if self.step.dice.randominit == "first" and cycle > 1:
                f.write("init = yesreadxyz\n")
                f.write(f"dens = {self.step.dice.dens:<8.4f}\n")
                f.write(f"vstep = {int(self.step.altsteps / 5)}\n")
            else:
                f.write("init = no\n")
                f.write(f"vstep = {int(self.step.dice.nstep[1] / 5)}\n")

            f.write("nstep = 5\n")
            f.write("mstop = 1\n")
            f.write("accum = no\n")
            f.write("iprint = 1\n")
            f.write("isave = 0\n")
            f.write("irdf = 0\n")

            seed = int(1e6 * random.random())
            f.write(f"seed = {seed}\n")

    def _make_npt_eq(self, proc_dir):
        file = Path(proc_dir, "NPT.eq")
        with open(file, "w") as f:
            f.write(f"title = {self.title} - NPT Production\n")
            f.write(f"ncores = {self.step.ncores}\n")
            f.write(f"ljname = {self.step.dice.ljname}\n")
            f.write(f"outname = {self.step.dice.outname}\n")

            mol_string = " ".join(str(x) for x in self.step.dice.nmol)
            f.write(f"nmol = {mol_string}\n")

            f.write(f"press = {self.step.dice.press}\n")
            f.write(f"temp = {self.step.dice.temp}\n")

            f.write(f"nstep = 5\n")

            f.write(f"vstep = {int(self.step.dice.nstep[2] / 5)}\n")
            f.write("init = no\n")
            f.write("mstop = 1\n")
            f.write("accum = no\n")
            f.write("iprint = 1\n")
            f.write(f"isave = {self.step.dice.isave}\n")
            f.write(f"irdf = {10 * self.step.nprocs}\n")

            seed = int(1e6 * random.random())
            f.write(f"seed = {seed}\n")

    def _make_potentials(self, proc_dir):
        fstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f}\n"

        file = Path(proc_dir, self.step.dice.ljname)
        with open(file, "w") as f:
            f.write(f"{self.step.dice.combrule}\n")
            f.write(f"{len(self.step.dice.nmol)}\n")

            nsites_qm = len(self.system.molecule[0].atom)
            f.write(f"{nsites_qm} {self.system.molecule[0].molname}\n")

            for atom in self.system.molecule[0].atom:
                f.write(
                    fstr.format(
                        atom.lbl,
                        atom.na,
                        atom.rx,
                        atom.ry,
                        atom.rz,
                        atom.chg,
                        atom.eps,
                        atom.sig,
                    )
                )

            for mol in self.system.molecule[1:]:
                f.write(f"{len(mol.atom)} {mol.molname}\n")
                for atom in mol.atom:
                    f.write(
                        fstr.format(
                            atom.lbl,
                            atom.na,
                            atom.rx,
                            atom.ry,
                            atom.rz,
                            atom.chg,
                            atom.eps,
                            atom.sig,
                        )
                    )

    def run_dice_file(self, cycle: int, proc: int, file_name: str):
        with open(Path(file_name), "r") as infile, open(
            Path(file_name + ".out"), "w"
        ) as outfile:
            if shutil.which("bash") is not None:
                exit_status = subprocess.call(
                    [
                        "bash",
                        "-c",
                        f"exec -a dice-step{cycle}-p{proc} {self.step.dice.progname} < {infile.name} > {outfile.name}",
                    ]
                )
            else:
                exit_status = subprocess.call(
                    self.step.dice.progname, stdin=infile, stdout=outfile
                )

        if exit_status != 0:
            raise RuntimeError(
                f"Dice process step{cycle:02d}-p{proc:02d} did not exit properly"
            )

        with open(Path(file_name + ".out"), "r") as outfile:
            flag = outfile.readlines()[DICE_FLAG_LINE].strip()
            if flag != DICE_END_FLAG:
                raise RuntimeError(
                    f"Dice process step{cycle:02d}-p{proc:02d} did not exit properly"
                )

        logger.info(f"Dice {file_name} - step{cycle:02d}-p{proc:02d} exited properly")
