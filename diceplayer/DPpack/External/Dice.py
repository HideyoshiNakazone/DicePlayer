import os
import shutil
import subprocess
import sys
from multiprocessing import Process, connection
from typing import Final, List, TextIO

import setproctitle
from diceplayer.DPpack.Utils.Misc import *
from diceplayer.DPpack.Utils.PTable import *
from diceplayer.DPpack.Utils.StepDTO import StepDTO
from diceplayer.DPpack.Utils.Validations import NotNull
from numpy import random

DICE_END_FLAG: Final[str] = "End of simulation"
DICE_FLAG_LINE: Final[int] = -2
UMAANG3_TO_GCM3: Final[float] = 1.6605

MAX_SEED: Final[int] = 4294967295


class Dice:

    title = "Diceplayer run"
    progname = "dice"

    nprocs: int = None
    randominit = "first"
    combrule = "*"

    temp = 300.0
    press = 1.0
    isave = 1000
    dens = None
    ljname = None
    outname = None
    nmol: List[int] = None
    nstep: List[int] = None
    upbuf = 360

    def __init__(self, infile: TextIO, outfile: TextIO) -> None:

        self.infile = infile
        self.outfile = outfile

    @NotNull(requiredArgs=["ncores", "nmol", "dens", "nstep", "ljname", "outname"])
    def updateKeywords(self, **data):
        self.__dict__.update(**data)

    def __new_density(self, cycle: int, proc: int) -> float:

        sim_dir = "simfiles"
        step_dir = "step{:02d}".format(cycle - 1)
        proc_dir = "p{:02d}".format(proc)
        path = sim_dir + os.sep + step_dir + os.sep + proc_dir
        file = path + os.sep + "last.xyz"

        if not os.path.isfile(file):
            sys.exit(
                "Error: cannot find the xyz file {} in main directory".format(file)
            )
        try:
            with open(file) as fh:
                xyzfile = fh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        box = xyzfile[1].split()
        volume = float(box[-3]) * float(box[-2]) * float(box[-1])

        total_mass = 0
        for i in range(len(self.step.molecule)):

            total_mass += self.step.molecule[i].total_mass * self.step.nmol[i]

        density = (total_mass / volume) * UMAANG3_TO_GCM3

        return density

    def __print_last_config(self, cycle: int, proc: int) -> None:

        sim_dir = "simfiles"
        step_dir = "step{:02d}".format(cycle)
        proc_dir = "p{:02d}".format(proc)
        path = sim_dir + os.sep + step_dir + os.sep + proc_dir
        file = path + os.sep + "phb.xyz"
        if not os.path.isfile(file):
            sys.exit("Error: cannot find the xyz file {}".format(file))
        try:
            with open(file) as fh:
                xyzfile = fh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        nsites = len(self.step.molecule[0].atom) * self.step.nmol[0]
        for i in range(1, len(self.step.nmol)):
            nsites += self.step.nmol[i] * len(self.step.molecule[i].atom)

        nsites += 2

        nsites *= -1
        xyzfile = xyzfile[nsites:]

        file = path + os.sep + "last.xyz"
        fh = open(file, "w")
        for line in xyzfile:
            fh.write(line)

    def __make_dice_inputs(self, cycle: int, proc: int) -> None:

        sim_dir = "simfiles"
        step_dir = "step{:02d}".format(cycle)
        proc_dir = "p{:02d}".format(proc)
        path = sim_dir + os.sep + step_dir + os.sep + proc_dir

        num = time.time()
        num = (num - int(num)) * 1e6

        num = int((num - int(num)) * 1e6)
        random.seed((os.getpid() * num) % (MAX_SEED + 1))

        if self.randominit == "first" and cycle > self.step.initcyc:
            last_step_dir = "step{:02d}".format(cycle - 1)
            last_path = sim_dir + os.sep + last_step_dir + os.sep + proc_dir
            xyzfile = last_path + os.sep + "last.xyz"
            self.__make_init_file(path, xyzfile)

        if len(self.nstep) == 2:

            self.__make_nvt_ter(cycle, path)
            self.__make_nvt_eq(path)

        elif len(self.nstep) == 3:

            if self.randominit == "first" and cycle > self.step.initcyc:
                self.dens = self.__new_density(cycle, proc)
            else:
                self.__make_nvt_ter(cycle, path)

            self.__make_npt_ter(cycle, path)
            self.__make_npt_eq(path)

        else:
            sys.exit("Error: bad number of entries for 'nstep'")

        self.__make_potential(path)

    def __make_nvt_ter(self, cycle: int, path: str) -> None:

        file = path + os.sep + "NVT.ter"
        try:
            fh = open(file, "w")
        except:
            sys.exit("Error: cannot open file {}".format(file))

        fh.write("title = {} - NVT Thermalization\n".format(self.title))
        fh.write("ncores = {}\n".format(self.ncores))
        fh.write("ljname = {}\n".format(self.ljname))
        fh.write("outname = {}\n".format(self.outname))

        string = " ".join(str(x) for x in self.nmol)
        fh.write("nmol = {}\n".format(string))

        fh.write("dens = {}\n".format(self.dens))
        fh.write("temp = {}\n".format(self.temp))

        if self.randominit == "first" and cycle > self.step.initcyc:
            fh.write("init = yesreadxyz\n")
            fh.write("nstep = {}\n".format(self.step.altsteps))
        else:
            fh.write("init = yes\n")
            fh.write("nstep = {}\n".format(self.nstep[0]))

        fh.write("vstep = 0\n")
        fh.write("mstop = 1\n")
        fh.write("accum = no\n")
        fh.write("iprint = 1\n")
        fh.write("isave = 0\n")
        fh.write("irdf = 0\n")

        seed = int(1e6 * random.random())
        fh.write("seed = {}\n".format(seed))
        fh.write("upbuf = {}".format(self.upbuf))

        fh.close()

    def __make_nvt_eq(self, path: str) -> None:

        file = path + os.sep + "NVT.eq"
        try:
            fh = open(file, "w")
        except:
            sys.exit("Error: cannot open file {}".format(file))

        fh.write("title = {} - NVT Production\n".format(self.title))
        fh.write("ncores = {}\n".format(self.ncores))
        fh.write("ljname = {}\n".format(self.ljname))
        fh.write("outname = {}\n".format(self.outname))

        string = " ".join(str(x) for x in self.nmol)
        fh.write("nmol = {}\n".format(string))

        fh.write("dens = {}\n".format(self.dens))
        fh.write("temp = {}\n".format(self.temp))
        fh.write("init = no\n")
        fh.write("nstep = {}\n".format(self.nstep[1]))
        fh.write("vstep = 0\n")
        fh.write("mstop = 1\n")
        fh.write("accum = no\n")
        fh.write("iprint = 1\n")
        fh.write("isave = {}\n".format(self.isave))
        fh.write("irdf = {}\n".format(10 * self.step.nprocs))

        seed = int(1e6 * random.random())
        fh.write("seed = {}\n".format(seed))

        fh.close()

    def __make_npt_ter(self, cycle: int, path: str) -> None:

        file = path + os.sep + "NPT.ter"
        try:
            fh = open(file, "w")
        except:
            sys.exit("Error: cannot open file {}".format(file))

        fh.write("title = {} - NPT Thermalization\n".format(self.title))
        fh.write("ncores = {}\n".format(self.ncores))
        fh.write("ljname = {}\n".format(self.ljname))
        fh.write("outname = {}\n".format(self.outname))

        string = " ".join(str(x) for x in self.nmol)
        fh.write("nmol = {}\n".format(string))

        fh.write("press = {}\n".format(self.press))
        fh.write("temp = {}\n".format(self.temp))

        if self.randominit == "first" and cycle > self.step.initcyc:
            fh.write("init = yesreadxyz\n")
            fh.write("dens = {:<8.4f}\n".format(self.dens))
            fh.write("vstep = {}\n".format(int(self.step.altsteps / 5)))
        else:
            fh.write("init = no\n")
            fh.write("vstep = {}\n".format(int(self.nstep[1] / 5)))

        fh.write("nstep = 5\n")
        fh.write("mstop = 1\n")
        fh.write("accum = no\n")
        fh.write("iprint = 1\n")
        fh.write("isave = 0\n")
        fh.write("irdf = 0\n")

        seed = int(1e6 * random.random())
        fh.write("seed = {}\n".format(seed))

        fh.close()

    def __make_npt_eq(self, path: str) -> None:

        file = path + os.sep + "NPT.eq"
        try:
            fh = open(file, "w")
        except:
            sys.exit("Error: cannot open file {}".format(file))

        fh.write("title = {} - NPT Production\n".format(self.title))
        fh.write("ncores = {}\n".format(self.ncores))
        fh.write("ljname = {}\n".format(self.ljname))
        fh.write("outname = {}\n".format(self.outname))

        string = " ".join(str(x) for x in self.nmol)
        fh.write("nmol = {}\n".format(string))

        fh.write("press = {}\n".format(self.press))
        fh.write("temp = {}\n".format(self.temp))

        fh.write("nstep = 5\n")

        fh.write("vstep = {}\n".format(int(self.nstep[2] / 5)))
        fh.write("init = no\n")
        fh.write("mstop = 1\n")
        fh.write("accum = no\n")
        fh.write("iprint = 1\n")
        fh.write("isave = {}\n".format(self.isave))
        fh.write("irdf = {}\n".format(10 * self.step.nprocs))

        seed = int(1e6 * random.random())
        fh.write("seed = {}\n".format(seed))

        fh.close()

    def __make_init_file(self, path: str, file: TextIO) -> None:

        if not os.path.isfile(file):
            sys.exit(
                "Error: cannot find the xyz file {} in main directory".format(file)
            )
        try:
            with open(file) as fh:
                xyzfile = fh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        nsites_mm = 0
        for i in range(1, len(self.step.nmol)):
            nsites_mm += self.step.nmol[i] * len(self.step.molecule[i].atom)

        nsites_mm *= -1

        xyzfile = xyzfile[nsites_mm:]

        file = path + os.sep + self.outname + ".xy"

        try:
            fh = open(file, "w", 1)
        except:
            sys.exit("Error: cannot open file {}".format(file))

        for atom in self.step.molecule[0].atom:
            fh.write(
                "{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(atom.rx, atom.ry, atom.rz)
            )

        for line in xyzfile:
            atom = line.split()
            rx = float(atom[1])
            ry = float(atom[2])
            rz = float(atom[3])
            fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(rx, ry, rz))

        fh.write("$end")

        fh.close()

    def __make_potential(self, path: str) -> None:

        fstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f}\n"

        file = path + os.sep + self.ljname
        try:
            fh = open(file, "w")
        except:
            sys.exit("Error: cannot open file {}".format(file))

        fh.write("{}\n".format(self.combrule))
        fh.write("{}\n".format(len(self.step.nmol)))

        nsites_qm = (
            len(self.step.molecule[0].atom)
            + len(self.step.molecule[0].ghost_atoms)
            + len(self.step.molecule[0].lp_atoms)
        )

        fh.write("{} {}\n".format(nsites_qm, self.step.molecule[0].molname))
        for atom in self.step.molecule[0].atom:
            fh.write(
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

        ghost_label = self.step.molecule[0].atom[-1].lbl + 1
        for i in self.step.molecule[0].ghost_atoms:
            fh.write(
                fstr.format(
                    ghost_label,
                    ghost_number,
                    self.step.molecule[0].atom[i].rx,
                    self.step.molecule[0].atom[i].ry,
                    self.step.molecule[0].atom[i].rz,
                    self.step.molecule[0].atom[i].chg,
                    0,
                    0,
                )
            )

        ghost_label += 1
        for lp in self.step.molecule[0].lp_atoms:
            fh.write(
                fstr.format(
                    ghost_label,
                    ghost_number,
                    lp["rx"],
                    lp["ry"],
                    lp["rz"],
                    lp["chg"],
                    0,
                    0,
                )
            )

        for mol in self.step.molecule[1:]:
            fh.write("{} {}\n".format(len(mol.atom), mol.molname))
            for atom in mol.atom:
                fh.write(
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

    def __make_proc_dir(self, cycle: int, proc: int) -> None:

        sim_dir = "simfiles"
        step_dir = "step{:02d}".format(cycle)
        proc_dir = "p{:02d}".format(proc)
        path = sim_dir + os.sep + step_dir + os.sep + proc_dir
        try:
            os.makedirs(path)
        except:
            sys.exit("Error: cannot make directory {}".format(path))

    def __run_dice(self, cycle: int, proc: int, fh: str) -> None:

        sim_dir = "simfiles"
        step_dir = "step{:02d}".format(cycle)
        proc_dir = "p{:02d}".format(proc)

        try:
            fh.write(
                "Simulation process {} initiated with pid {}\n".format(
                    sim_dir + os.sep + step_dir + os.sep + proc_dir, os.getpid()
                )
            )

        except Exception as err:
            print("I/O error({0}): {1}".format(err))

        path = sim_dir + os.sep + step_dir + os.sep + proc_dir
        working_dir = os.getcwd()
        os.chdir(path)

        if len(self.nstep) == 2:

            if self.randominit == "first" and cycle > self.step.initcyc:
                string_tmp = "previous"
            else:
                string_tmp = "random"

            string = "(from " + string_tmp + " configuration)"
            fh.write(
                "p{:02d}> NVT thermalization finished {} on {}\n".format(
                    proc, string, date_time()
                )
            )

            infh = open("NVT.ter")
            outfh = open("NVT.ter.out", "w")

            if shutil.which("bash") != None:
                exit_status = subprocess.call(
                    [
                        "bash",
                        "-c",
                        "exec -a dice-step{}-p{} {} < {} > {}".format(
                            cycle, proc, self.progname, infh.name, outfh.name
                        ),
                    ]
                )
            else:
                exit_status = subprocess.call(
                    self.progname, stin=infh.name, stout=outfh.name
                )

            infh.close()
            outfh.close()

            if os.getppid() == 1:
                sys.exit()

            if exit_status != 0:
                sys.exit(
                    "Dice process step{:02d}-p{:02d} did not exit properly".format(
                        cycle, proc
                    )
                )
            else:
                outfh = open("NVT.ter.out")
                flag = outfh.readlines()[DICE_FLAG_LINE].strip()
                outfh.close()
                if flag != DICE_END_FLAG:
                    sys.exit(
                        "Dice process step{:02d}-p{:02d} did not exit properly".format(
                            cycle, proc
                        )
                    )

            fh.write(
                "p{:02d}> NVT production initiated on {}\n".format(proc, date_time())
            )

            infh = open("NVT.eq")
            outfh = open("NVT.eq.out", "w")

            if shutil.which("bash") != None:
                exit_status = subprocess.call(
                    [
                        "bash",
                        "-c",
                        "exec -a dice-step{}-p{} {} < {} > {}".format(
                            cycle, proc, self.progname, infh.name, outfh.name
                        ),
                    ]
                )
            else:
                exit_status = subprocess.call(
                    self.progname, stin=infh.name, stout=outfh.name
                )

            infh.close()
            outfh.close()

            if os.getppid() == 1:
                sys.exit()

            if exit_status != 0:
                sys.exit(
                    "Dice process step{:02d}-p{:02d} did not exit properly".format(
                        cycle, proc
                    )
                )
            else:
                outfh = open("NVT.eq.out")
                flag = outfh.readlines()[DICE_FLAG_LINE].strip()
                outfh.close()
                if flag != DICE_END_FLAG:
                    sys.exit(
                        "Dice process step{:02d}-p{:02d} did not exit properly".format(
                            cycle, proc
                        )
                    )

            fh.write(
                "p{:02d}> ----- NVT production finished on {}\n".format(
                    proc, date_time()
                )
            )

        elif len(self.nstep) == 3:
            if (
                self.randominit == "always"
                or (self.randominit == "first" and cycle == 1)
                or self.continued
            ):
                string = "(from random configuration)"
                fh.write(
                    "p{:02d}> NVT thermalization initiated {} on {}\n".format(
                        proc, string, date_time()
                    )
                )
                infh = open("NVT.ter")
                outfh = open("NVT.ter.out", "w")

                if shutil.which("bash") != None:
                    exit_status = subprocess.call(
                        [
                            "bash",
                            "-c",
                            "exec -a dice-step{}-p{} {} < {} > {}".format(
                                cycle, proc, self.progname, infh.name, outfh.name
                            ),
                        ]
                    )
                else:
                    exit_status = subprocess.call(
                        self.progname, stin=infh.name, stout=outfh.name
                    )

                infh.close()
                outfh.close()

                if os.getppid() == 1:
                    sys.exit()

                if exit_status != 0:
                    sys.exit(
                        "Dice process step{:02d}-p{:02d} did not exit properly".format(
                            cycle, proc
                        )
                    )
                else:
                    outfh = open("NVT.ter.out")
                    flag = outfh.readlines()[DICE_FLAG_LINE].strip()
                    outfh.close()
                    if flag != DICE_END_FLAG:
                        sys.exit(
                            "Dice process step{:02d}-p{:02d} did not exit properly".format(
                                cycle, proc
                            )
                        )

            if not self.randominit == "always" or (
                (self.randominit == "first" and cycle > self.step.initcyc)
            ):
                string = " (from previous configuration) "
            else:
                string = " "
            fh.write(
                "p{:02d}> NPT thermalization finished {} on {}\n".format(
                    proc, string, date_time()
                )
            )

            infh = open("NPT.ter")
            outfh = open("NPT.ter.out", "w")

            if shutil.which("bash") != None:
                exit_status = subprocess.call(
                    [
                        "bash",
                        "-c",
                        "exec -a dice-step{}-p{} {} < {} > {}".format(
                            cycle, proc, self.progname, infh.name, outfh.name
                        ),
                    ]
                )
            else:
                exit_status = subprocess.call(
                    self.progname, stin=infh.name, stout=outfh.name
                )

            infh.close()
            outfh.close()

            if os.getppid() == 1:
                sys.exit()

            if exit_status != 0:
                sys.exit(
                    "Dice process step{:02d}-p{:02d} did not exit properly".format(
                        cycle, proc
                    )
                )
            else:
                outfh = open("NPT.ter.out")
                flag = outfh.readlines()[DICE_FLAG_LINE].strip()
                outfh.close()
                if flag != DICE_END_FLAG:
                    sys.exit(
                        "Dice process step{:02d}-p{:02d} did not exit properly".format(
                            cycle, proc
                        )
                    )

            fh.write(
                "p{:02d}> NPT production initiated on {}\n".format(proc, date_time())
            )

            infh = open("NPT.eq")
            outfh = open("NPT.eq.out", "w")

            if shutil.which("bash") != None:
                exit_status = subprocess.call(
                    [
                        "bash",
                        "-c",
                        "exec -a dice-step{}-p{} {} < {} > {}".format(
                            cycle, proc, self.progname, infh.name, outfh.name
                        ),
                    ]
                )
            else:
                exit_status = subprocess.call(
                    self.progname, stin=infh.name, stout=outfh.name
                )

            infh.close()
            outfh.close()

            if os.getppid() == 1:
                sys.exit()

            if exit_status != 0:
                sys.exit(
                    "Dice process step{:02d}-p{:02d} did not exit properly".format(
                        cycle, proc
                    )
                )
            else:
                outfh = open("NPT.eq.out")
                flag = outfh.readlines()[DICE_FLAG_LINE].strip()
                outfh.close()
                if flag != DICE_END_FLAG:
                    sys.exit(
                        "Dice process step{:02d}-p{:02d} did not exit properly".format(
                            cycle, proc
                        )
                    )

            fh.write(
                "p{:02d}> ----- NPT production finished on {}\n".format(
                    proc, date_time()
                )
            )

        os.chdir(working_dir)

    def __simulation_process(self, cycle: int, proc: int):
        setproctitle.setproctitle("diceplayer-step{:0d}-p{:0d}".format(cycle, proc))

        try:
            self.__make_proc_dir(cycle, proc)
            self.__make_dice_inputs(cycle, proc)
            self.__run_dice(cycle, proc, self.outfile)
        except Exception as err:
            sys.exit(err)

    def configure(self, step: StepDTO):
        self.step = step

    def start(self, cycle: int) -> None:

        procs = []
        sentinels = []

        for proc in range(1, self.step.nprocs + 1):

            p = Process(target=self.__simulation_process, args=(cycle, proc))
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

        for proc in range(1, self.step.nprocs + 1):
            self.__print_last_config(cycle, proc)

    def reset(self):
        del self.step
