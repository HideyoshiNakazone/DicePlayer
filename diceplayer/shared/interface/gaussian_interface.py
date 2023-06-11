from __future__ import annotations

from diceplayer import logger
from diceplayer.shared.config.player_config import PlayerConfig
from diceplayer.shared.environment.atom import Atom
from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.environment.system import System
from diceplayer.shared.interface import Interface
from diceplayer.shared.utils.misc import date_time
from diceplayer.shared.utils.ptable import atomsymb

import numpy as np
from nptyping import NDArray

import os
import shutil
import subprocess
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Tuple


class GaussianInterface(Interface):
    def configure(self, step_dto: PlayerConfig, system: System):
        self.system = system
        self.step = step_dto

    def start(self, cycle: int) -> Dict[str, NDArray]:
        self._make_qm_dir(cycle)

        if cycle > 1:
            self._copy_chk_file_from_previous_step(cycle)

        asec_charges = self.populate_asec_vdw(cycle)
        self._make_gaussian_input_file(cycle, asec_charges)

        self._run_gaussian(cycle)
        self._run_formchk(cycle)

        return_value = {}
        if self.step.opt:
            # return_value['position'] = np.array(
            #     self._run_optimization(cycle)
            # )
            raise NotImplementedError("Optimization not implemented yet.")

        else:
            return_value["charges"] = np.array(self._read_charges_from_fchk(cycle))

        return return_value

    def reset(self):
        del self.step
        del self.system

    def _make_qm_dir(self, cycle: int):
        qm_dir_path = Path(self.step.simulation_dir, f"step{cycle:02d}", "qm")
        if not qm_dir_path.exists():
            qm_dir_path.mkdir()

    def _copy_chk_file_from_previous_step(self, cycle: int):
        current_chk_file_path = Path(
            self.step.simulation_dir, f"step{cycle:02d}", "qm", f"asec.chk"
        )
        if current_chk_file_path.exists():
            raise FileExistsError(f"File {current_chk_file_path} already exists.")

        previous_chk_file_path = Path(
            self.step.simulation_dir, f"step{(cycle - 1):02d}", "qm", f"asec.chk"
        )
        if not previous_chk_file_path.exists():
            raise FileNotFoundError(f"File {previous_chk_file_path} does not exist.")

        shutil.copy(previous_chk_file_path, current_chk_file_path)

    def populate_asec_vdw(self, cycle: int) -> list[dict]:
        norm_factor = self._calculate_norm_factor()

        nsitesref = len(self.system.molecule[0].atom)

        nsites_total = self._calculate_total_number_of_sites(nsitesref)

        proc_charges = []
        for proc in range(1, self.step.nprocs + 1):
            proc_charges.append(self._read_charges_from_last_step(cycle, proc))

        asec_charges, thickness, picked_mols = self._evaluate_proc_charges(
            nsites_total, proc_charges
        )

        logger.info(
            f"In average, {(sum(picked_mols) / norm_factor):^7.2f} molecules\n"
            f"were selected from each of the {len(picked_mols)} configurations\n"
            f"of the production simulations to form the ASEC, comprising a shell with\n"
            f"minimum thickness of {(sum(thickness) / norm_factor):>6.2f} Angstrom\n"
        )

        for charge in asec_charges:
            charge["chg"] = charge["chg"] / norm_factor

        return asec_charges

    def _calculate_norm_factor(self) -> int:
        if self.step.dice.nstep[-1] % self.step.dice.isave == 0:
            nconfigs = round(self.step.dice.nstep[-1] / self.step.dice.isave)
        else:
            nconfigs = int(self.step.dice.nstep[-1] / self.step.dice.isave)

        return nconfigs * self.step.nprocs

    def _calculate_total_number_of_sites(self, nsitesref) -> int:
        nsites_total = self.step.dice.nmol[0] * nsitesref
        for i in range(1, len(self.step.dice.nmol)):
            nsites_total += self.step.dice.nmol[i] * len(self.system.molecule[i].atom)

        return nsites_total

    def _read_charges_from_last_step(self, cycle: int, proc: int) -> list[str]:
        last_xyz_file_path = Path(
            self.step.simulation_dir, f"step{cycle:02d}", f"p{proc:02d}", "last.xyz"
        )
        if not last_xyz_file_path.exists():
            raise FileNotFoundError(f"File {last_xyz_file_path} does not exist.")

        with open(last_xyz_file_path, "r") as last_xyz_file:
            lines = last_xyz_file.readlines()

        return lines

    def _evaluate_proc_charges(
        self, total_nsites: int, proc_charges: list[list[str]]
    ) -> Tuple[List[Dict[str, float | Any]], List[float], List[int]]:
        asec_charges = []

        thickness = []
        picked_mols = []

        for charges in proc_charges:
            charges_nsites = int(charges.pop(0))
            if int(charges_nsites) != total_nsites:
                raise ValueError(
                    f"Number of sites does not match total number of sites."
                )

            thickness.append(self._calculate_proc_thickness(charges))
            nsites_ref_mol = len(self.system.molecule[0].atom)
            charges = charges[nsites_ref_mol:]

            mol_count = 0
            for type in range(len(self.step.dice.nmol)):
                if type == 0:
                    # Reference Molecule must be ignored from type 0
                    nmols = self.step.dice.nmol[type] - 1
                else:
                    nmols = self.step.dice.nmol[type]

                for mol in range(nmols):
                    new_molecule = Molecule("ASEC TMP MOLECULE")
                    for site in range(len(self.system.molecule[type].atom)):
                        line = charges.pop(0).split()

                        if (
                            line[0].title()
                            != atomsymb[
                                self.system.molecule[type].atom[site].na
                            ].strip()
                        ):
                            raise SyntaxError(
                                f"Error: Invalid Dice Output. Atom type does not match."
                            )

                        new_molecule.add_atom(
                            Atom(
                                self.system.molecule[type].atom[site].lbl,
                                self.system.molecule[type].atom[site].na,
                                float(line[1]),
                                float(line[2]),
                                float(line[3]),
                                self.system.molecule[type].atom[site].chg,
                                self.system.molecule[type].atom[site].eps,
                                self.system.molecule[type].atom[site].sig,
                            )
                        )

                    distance = self.system.molecule[0].minimum_distance(new_molecule)

                    if distance < thickness[-1]:
                        for atom in new_molecule.atom:
                            asec_charges.append(
                                {
                                    "lbl": atomsymb[atom.na],
                                    "rx": atom.rx,
                                    "ry": atom.ry,
                                    "rz": atom.rz,
                                    "chg": atom.chg,
                                }
                            )
                        mol_count += 1

            picked_mols.append(mol_count)

        return asec_charges, thickness, picked_mols

    def _calculate_proc_thickness(self, charges: list[str]) -> float:
        box = charges.pop(0).split()[-3:]
        box = [float(box[0]), float(box[1]), float(box[2])]
        sizes = self.system.molecule[0].sizes_of_molecule()

        return min(
            [
                (box[0] - sizes[0]) / 2,
                (box[1] - sizes[1]) / 2,
                (box[2] - sizes[2]) / 2,
            ]
        )

    def _make_gaussian_input_file(self, cycle: int, asec_charges: list[dict]) -> None:
        gaussian_input_file_path = Path(
            self.step.simulation_dir, f"step{cycle:02d}", "qm", f"asec.gjf"
        )

        with open(gaussian_input_file_path, "w") as gaussian_input_file:
            gaussian_input_file.writelines(
                self._generate_gaussian_input(cycle, asec_charges)
            )

    def _generate_gaussian_input(
        self, cycle: int, asec_charges: list[dict]
    ) -> list[str]:
        gaussian_input = ["%Chk=asec.chk\n"]

        if self.step.mem is not None:
            gaussian_input.append(f"%Mem={self.step.mem}GB\n")

        gaussian_input.append(f"%Nprocs={self.step.nprocs * self.step.ncores}\n")

        kwords_line = f"#P {self.step.gaussian.level}"

        if self.step.gaussian.keywords:
            kwords_line += " " + self.step.gaussian.keywords

        if self.step.opt == "yes":
            kwords_line += " Force"

        kwords_line += " NoSymm"
        kwords_line += f" Pop={self.step.gaussian.pop} Density=Current"

        if cycle > 1:
            kwords_line += " Guess=Read"

        gaussian_input.append(textwrap.fill(kwords_line, 90))
        gaussian_input.append("\n")

        gaussian_input.append("\nForce calculation - Cycle number {}\n".format(cycle))
        gaussian_input.append("\n")
        gaussian_input.append(
            f"{self.step.gaussian.chgmult[0]},{self.step.gaussian.chgmult[1]}\n"
        )

        for atom in self.system.molecule[0].atom:
            symbol = atomsymb[atom.na]
            gaussian_input.append(
                "{:<2s}    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
                    symbol, atom.rx, atom.ry, atom.rz
                )
            )

        gaussian_input.append("\n")

        for charge in asec_charges:
            gaussian_input.append(
                "{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
                    charge["rx"], charge["ry"], charge["rz"], charge["chg"]
                )
            )

        gaussian_input.append("\n")

        return gaussian_input

    def _run_gaussian(self, cycle: int) -> None:
        qm_dir = Path(self.step.simulation_dir, f"step{(cycle):02d}", "qm")

        working_dir = os.getcwd()
        os.chdir(qm_dir)

        infile = "asec.gjf"

        operation = None
        if self.step.opt:
            operation = "forces"
        else:
            operation = "charges"

        logger.info(
            f"Calculation of {operation} initiated with Gaussian on {date_time()}\n"
        )

        if shutil.which("bash") is not None:
            exit_status = subprocess.call(
                [
                    "bash",
                    "-c",
                    "exec -a {}-step{} {} {}".format(
                        self.step.gaussian.qmprog,
                        cycle,
                        self.step.gaussian.qmprog,
                        infile,
                    ),
                ]
            )
        else:
            exit_status = subprocess.call([self.step.gaussian.qmprog, infile])

        if exit_status != 0:
            raise SystemError("Gaussian process did not exit properly")

        logger.info(f"Calculation of {operation} finished on {date_time()}")

        os.chdir(working_dir)

    def _run_formchk(self, cycle: int):
        qm_dir = Path(self.step.simulation_dir, f"step{(cycle):02d}", "qm")

        work_dir = os.getcwd()
        os.chdir(qm_dir)

        logger.info("Formatting the checkpoint file... \n")

        exit_status = subprocess.call(
            ["formchk", "asec.chk"], stdout=subprocess.DEVNULL
        )

        if exit_status != 0:
            raise SystemError("Formchk process did not exit properly")

        logger.info("Done\n")

        os.chdir(work_dir)

    def _read_charges_from_fchk(self, cycle: int):
        fchk_file_path = Path("simfiles", f"step{cycle:02d}", "qm", "asec.fchk")
        with open(fchk_file_path) as fchk:
            fchkfile = fchk.readlines()

        if self.step.gaussian.pop in ["chelpg", "mk"]:
            CHARGE_FLAG = "ESP Charges"
        else:
            CHARGE_FLAG = "ESP Charges"

        start = fchkfile.pop(0).strip()
        while start.find(CHARGE_FLAG) != 0:  # expression in begining of line
            start = fchkfile.pop(0).strip()

        charges: List[float] = []
        while len(charges) < len(self.system.molecule[0].atom):
            charges.extend([float(x) for x in fchkfile.pop(0).split()])

        return charges
