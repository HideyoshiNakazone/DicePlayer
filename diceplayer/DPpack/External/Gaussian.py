import os
import shutil
import subprocess
import sys
import textwrap
from typing import TextIO

import numpy as np

from diceplayer.DPpack.Environment.Atom import Atom
from diceplayer.DPpack.Environment.Molecule import Molecule
from diceplayer.DPpack.Utils.Misc import *
from diceplayer.DPpack.Utils.PTable import *
from diceplayer.DPpack.Utils.StepDTO import StepDTO
from diceplayer.DPpack.Utils.Validations import NotNull


class Gaussian:

    qmprog = "g09"
    mem = None
    keywords = None
    chgmult = [0, 1]
    gmiddle = None  # In each case, if a filename is given, its content will be
    gbottom = None  # inserted in the gaussian input
    pop = "chelpg"
    level = None

    def __init__(self) -> None:
        pass

    @NotNull(requiredArgs=["level"])
    def updateKeywords(self, **data):
        self.__dict__.update(data)

    def run_formchk(self, cycle: int, fh: TextIO):

        simdir = "simfiles"
        stepdir = "step{:02d}".format(cycle)
        path = simdir + os.sep + stepdir + os.sep + "qm"

        work_dir = os.getcwd()
        os.chdir(path)

        fh.write("Formatting the checkpoint file... \n")

        exit_status = subprocess.call(["formchk", "asec.chk"], stdout=fh)

        fh.write("Done\n")

        os.chdir(work_dir)

    def read_forces_fchk(self, file: str, fh: TextIO) -> np.ndarray:

        forces = []
        try:
            with open(file) as tmpfh:
                fchkfile = tmpfh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        start = fchkfile.pop(0).strip()
        while start.find("Cartesian Gradient") != 0:  # expression in begining of line
            start = fchkfile.pop(0).strip()

        degrees = 3 * len(self.step.molecule[0].atom)
        count = 0
        while len(forces) < degrees:
            values = fchkfile.pop(0).split()
            forces.extend([float(x) for x in values])
            count += len(values)
            if count >= degrees:
                forces = forces[:degrees]
                break

        gradient = np.array(forces)

        fh.write("\nGradient read from file {}:\n".format(file))
        fh.write(
            "-----------------------------------------------------------------------\n"
            "Center     Atomic                     Forces (Hartree/Bohr)\n"
            "Number     Number              X                Y                Z\n"
            "-----------------------------------------------------------------------\n"
        )
        for i in range(len(self.step.molecule[0].atom)):
            fh.write(
                "  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
                    i + 1,
                    self.step.molecule[0].atom[i].na,
                    forces.pop(0),
                    forces.pop(0),
                    forces.pop(0),
                )
            )

        fh.write(
            "-----------------------------------------------------------------------\n"
        )

        force_max = np.amax(np.absolute(gradient))
        force_rms = np.sqrt(np.mean(np.square(gradient)))

        fh.write(
            "  Max Force = {:>14.9f}      RMS Force = {:>14.9f}\n\n".format(
                force_max, force_rms
            )
        )

        return gradient

    def read_hessian_fchk(self, file: str) -> np.ndarray:

        force_const = []
        try:
            with open(file) as tmpfh:
                fchkfile = tmpfh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        start = fchkfile.pop(0).strip()
        while start.find("Cartesian Force Constants") != 0:
            start = fchkfile.pop(0).strip()

        degrees = 3 * len(self.step.molecule[0].atom)
        last = round(degrees * (degrees + 1) / 2)
        count = 0

        while len(force_const) < last:

            value = fchkfile.pop(0).split()
            force_const.extend([float(x) for x in value])

        # while len(force_const) < last:

        # 	values = fchkfile.pop(0).split()
        # 	force_const.extend([ float(x) for x in values ])
        # 	count += len(values)
        # 	if count >= last:
        # 		force_const = force_const[:last]
        # 		break

        hessian = np.zeros((degrees, degrees))
        for i in range(degrees):
            for j in range(i + 1):
                hessian[i, j] = force_const.pop(0)
                hessian[j, i] = hessian[i, j]

        return hessian

    def read_hessian_log(self, file: str) -> np.ndarray:

        try:
            with open(file) as tmpfh:
                logfile = tmpfh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        start = logfile.pop(0).strip()
        while start.find("The second derivative matrix:") != 0:
            start = logfile.pop(0).strip()

        degrees = 3 * len(self.step.molecule[0].atom)
        hessian = np.zeros((degrees, degrees))

        k = 0
        while k < degrees:
            logfile.pop(0)
            for i in range(k, degrees):
                values = logfile.pop(0).split()[1:]
                for j in range(k, min(i + 1, k + 5)):
                    hessian[i, j] = float(values.pop(0))
                    hessian[j, i] = hessian[i, j]
            k += 5

        return hessian

    def print_grad_hessian(
        self, cycle: int, cur_gradient: np.ndarray, hessian: np.ndarray
    ) -> None:

        try:
            fh = open("grad_hessian.dat", "w")
        except:
            sys.exit("Error: cannot open file grad_hessian.dat")

        fh.write("Optimization cycle: {}\n".format(cycle))
        fh.write("Cartesian Gradient\n")
        degrees = 3 * len(self.step.molecule[0].atom)
        for i in range(degrees):
            fh.write(" {:>11.8g}".format(cur_gradient[i]))
            if (i + 1) % 5 == 0 or i == degrees - 1:
                fh.write("\n")

        fh.write("Cartesian Force Constants\n")
        n = int(np.sqrt(2 * degrees))
        last = degrees * (degrees + 1) / 2
        count = 0
        for i in range(n):
            for j in range(i + 1):
                count += 1
                fh.write(" {:>11.8g}".format(hessian[i, j]))
                if count % 5 == 0 or count == last:
                    fh.write("\n")

        fh.close()

    # Change the name to make_gaussian_input
    def make_gaussian_input(self, cycle: int, asec_charges=None) -> None:

        simdir = "simfiles"
        stepdir = "step{:02d}".format(cycle)
        path = simdir + os.sep + stepdir + os.sep + "qm"

        file = path + os.sep + "asec.gjf"

        try:
            fh = open(file, "w")
        except:
            sys.exit("Error: cannot open file {}".format(file))

        fh.write("%Chk=asec.chk\n")
        if self.mem != None:
            fh.write("%Mem={}MB\n".format(self.mem))
        fh.write("%Nprocs={}\n".format(self.step.nprocs * self.step.ncores))

        kword_line = "#P " + str(self.level)

        if self.keywords != None:
            kword_line += " " + self.keywords

        if self.step.opt == "yes":
            kword_line += " Force"

        # kword_line += " Charge"
        kword_line += " NoSymm"
        kword_line += " Pop={} Density=Current".format(self.pop)

        if cycle > 1:
            kword_line += " Guess=Read"

        fh.write(textwrap.fill(kword_line, 90))
        fh.write("\n")

        fh.write("\nForce calculation - Cycle number {}\n".format(cycle))
        fh.write("\n")
        fh.write("{},{}\n".format(self.chgmult[0], self.chgmult[1]))

        for atom in self.step.molecule[0].atom:
            symbol = atomsymb[atom.na]
            fh.write(
                "{:<2s}    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
                    symbol, atom.rx, atom.ry, atom.rz
                )
            )

        # ## If also performing charge fit in the same calculation
        # if cycle >= self.player.switchcyc:
        # 	for ghost in ghost_atoms:
        # 		fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
        # 											ghost['rx'], ghost['ry'], ghost['rz']))

        # 	for lp in lp_atoms:
        # 		fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
        # 														lp['rx'], lp['ry'], lp['rz']))

        # fh.write("\n")

        # If gmiddle file was informed, write its contents in asec.gjf
        # if self.gmiddle != None:
        # 	if not os.path.isfile(self.gmiddle):
        # 		sys.exit("Error: cannot find file {} in main directory".format(
        # 																self.gmiddle))
        # 	try:
        # 		with open(self.gmiddle) as gmiddlefile:
        # 			gmiddle = gmiddlefile.readlines()
        # 	except:
        # 		sys.exit("Error: cannot open file {}".format(self.gmiddle))

        # 	for line in gmiddle:
        # 		fh.write(line)

        # 	fh.write("\n")

        # ## Write the ASEC:
        # for charge in asec_charges:
        # 	fh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
        # 							charge['rx'], charge['ry'], charge['rz'], charge['chg']))

        fh.write("\n")

        # ## If gbottom file was informed, write its contents in asec.gjf
        # if self.gbottom != None:
        # 	if not os.path.isfile(self.gbottom):
        # 		sys.exit("Error: cannot find file {} in main directory".format(
        # 																self.gbottom))
        # 	try:
        # 		with open(self.gbottom) as gbottomfile:
        # 			gbottom = gbottomfile.readlines()
        # 	except:
        # 		sys.exit("Error: cannot open file {}".format(self.gbottom))

        # 	for line in gbottom:
        # 		fh.write(line)

        # fh.write("\n")

        # fh.close()

    def read_charges(self, file: str, fh: TextIO) -> None:

        try:
            with open(file) as tmpfh:
                glogfile = tmpfh.readlines()
        except:
            sys.exit("Error: cannot open file {}".format(file))

        start = glogfile.pop(0).strip()
        while start != "Fitting point charges to electrostatic potential":
            start = glogfile.pop(0).strip()

        glogfile = glogfile[3:]  # Consume 3 more lines

        fh.write("\nAtomic charges:\n")
        fh.write("------------------------------------\n")
        for atom in self.step.molecule[0].atom:
            line = glogfile.pop(0).split()
            atom_str = line[1]
            charge = float(line[2])
            atom.chg = charge
            fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))

        # if self.pop == "chelpg":
        # 	for ghost in ghost_atoms:
        # 		line = glogfile.pop(0).split()
        # 		atom_str = line[1]
        # 		charge = float(line[2])
        # 		ghost['chg'] = charge
        # 		fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))

        # 	for lp in lp_atoms:
        # 		line = glogfile.pop(0).split()
        # 		atom_str = line[1]
        # 		charge = float(line[2])
        # 		lp['chg'] = charge
        # 		fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))

        fh.write("------------------------------------\n")

    def run_gaussian(self, cycle: int, type: str, fh: TextIO) -> None:

        simdir = "simfiles"
        stepdir = "step{:02d}".format(cycle)
        path = simdir + os.sep + stepdir + os.sep + "qm"
        work_dir = os.getcwd()
        os.chdir(path)

        # if type == "force":
        # 	infile = "asec.gjf"
        # elif type == "charge":
        # 	infile = "asec2.gjf"

        infile = "asec.gjf"

        fh.write(
            "\nCalculation of {}s initiated with Gaussian on {}\n".format(
                type, date_time()
            )
        )

        if shutil.which("bash") != None:
            exit_status = subprocess.call(
                [
                    "bash",
                    "-c",
                    "exec -a {}-step{} {} {}".format(
                        self.qmprog, cycle, self.qmprog, infile
                    ),
                ]
            )
        else:
            exit_status = subprocess.call([self.qmprog, infile])

        if exit_status != 0:
            sys.exit("Gaussian process did not exit properly")

        fh.write("Calculation of {}s finished on {}\n".format(type, date_time()))

        os.chdir(work_dir)

    # def calculate_step(
    #     self, cycle: int, gradient: np.ndarray, hessian: np.ndarray
    # ) -> np.ndarray:

    #     invhessian = np.linalg.inv(hessian)
    #     pre_step = -1 * np.matmul(invhessian, gradient.T).T
    #     maxstep = np.amax(np.absolute(pre_step))
    #     factor = min(1, self.player.maxstep / maxstep)
    #     step = factor * pre_step

    #     self.outfile.write("\nCalculated step-{}:\n".format(cycle))
    #     pre_step_list = pre_step.tolist()

    #     self.outfile.write(
    #         "-----------------------------------------------------------------------\n"
    #         "Center     Atomic                          Step (Bohr)\n"
    #         "Number     Number              X                Y                Z\n"
    #         "-----------------------------------------------------------------------\n"
    #     )
    #     for i in range(len(self.system.molecule[0].atom)):
    #         self.outfile.write(
    #             "  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
    #                 i + 1,
    #                 self.system.molecule[0].atom[i].na,
    #                 pre_step_list.pop(0),
    #                 pre_step_list.pop(0),
    #                 pre_step_list.pop(0),
    #             )
    #         )

    #     self.outfile.write(
    #         "-----------------------------------------------------------------------\n"
    #     )

    #     self.outfile.write("Maximum step is {:>11.6}\n".format(maxstep))
    #     self.outfile.write("Scaling factor = {:>6.4f}\n".format(factor))
    #     self.outfile.write("\nFinal step (Bohr):\n")
    #     step_list = step.tolist()

    #     self.outfile.write(
    #         "-----------------------------------------------------------------------\n"
    #         "Center     Atomic                          Step (Bohr)\n"
    #         "Number     Number              X                Y                Z\n"
    #         "-----------------------------------------------------------------------\n"
    #     )
    #     for i in range(len(self.system.molecule[0].atom)):
    #         self.outfile.write(
    #             "  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
    #                 i + 1,
    #                 self.system.molecule[0].atom[i].na,
    #                 step_list.pop(0),
    #                 step_list.pop(0),
    #                 step_list.pop(0),
    #             )
    #         )

    #     self.outfile.write(
    #         "-----------------------------------------------------------------------\n"
    #     )

    #     step_max = np.amax(np.absolute(step))
    #     step_rms = np.sqrt(np.mean(np.square(step)))

    #     self.outfile.write(
    #         "  Max Step = {:>14.9f}      RMS Step = {:>14.9f}\n\n".format(
    #             step_max, step_rms
    #         )
    #     )

    #     return step

    def configure(self, step: StepDTO):

        self.step = step

    def start(self, cycle: int, outfile: TextIO, readhessian: str) -> np.ndarray:

        make_qm_dir(cycle)

        if self.qmprog in ("g03", "g09", "g16"):

            if cycle > 1:

                src = (
                    "simfiles"
                    + os.sep
                    + "step{:02d}".format(cycle - 1)
                    + os.sep
                    + "qm"
                    + os.sep
                    + "asec.chk"
                )
                dst = (
                    "simfiles"
                    + os.sep
                    + "step{:02d}".format(cycle)
                    + os.sep
                    + "qm"
                    + os.sep
                    + "asec.chk"
                )
                shutil.copyfile(src, dst)

            self.make_gaussian_input(cycle)
            self.run_gaussian(cycle, "force", outfile)
            self.run_formchk(cycle, outfile)

            ## Read the gradient
            file = (
                "simfiles"
                + os.sep
                + "step{:02d}".format(cycle)
                + os.sep
                + "qm"
                + os.sep
                + "asec.fchk"
            )

            try:
                gradient
                old_gradient = gradient
            except:
                pass

            gradient = self.read_forces_fchk(file, outfile)

            # If 1st step, read the hessian
            if cycle == 1:

                if readhessian == "yes":

                    file = "grad_hessian.dat"
                    outfile.write(
                        "\nReading the hessian matrix from file {}\n".format(file)
                    )
                    hessian = self.read_hessian_log(file)

                else:

                    file = (
                        "simfiles"
                        + os.sep
                        + "step01"
                        + os.sep
                        + "qm"
                        + os.sep
                        + "asec.fchk"
                    )
                    outfile.write(
                        "\nReading the hessian matrix from file {}\n".format(file)
                    )
                    hessian = self.read_hessian_fchk(file)

            # From 2nd step on, update the hessian
            else:
                outfile.write("\nUpdating the hessian matrix using the BFGS method... ")
                hessian = self.step.molecule[0].update_hessian(
                    step, gradient, old_gradient, hessian
                )
                outfile.write("Done\n")

            # Save gradient and hessian
            self.print_grad_hessian(cycle, gradient, hessian)

            # Calculate the step and update the position
            step = self.calculate_step(cycle, gradient, hessian)

            position += step

            ## If needed, calculate the charges
            if cycle < self.step.switchcyc:

                # internal.gaussian.make_charge_input(cycle, asec_charges)
                self.run_gaussian(cycle, "charge", outfile)

                file = (
                    "simfiles"
                    + os.sep
                    + "step{:02d}".format(cycle)
                    + os.sep
                    + "qm"
                    + os.sep
                    + "asec2.log"
                )
                self.read_charges(file, outfile)
            else:
                file = (
                    "simfiles"
                    + os.sep
                    + "step{:02d}".format(cycle)
                    + os.sep
                    + "qm"
                    + os.sep
                    + "asec.log"
                )
                self.read_charges(file, outfile)

            ## Print new info for molecule[0]
            self.outfile.write("\nNew values for molecule type 1:\n\n")
            self.step.molecule[0].print_mol_info(outfile)

            ##
            ##  Molcas block
            ##
            # if player['qmprog'] == "molcas":

        # elif player['opt'] == "ts":

        ##
        ##  Gaussian block
        ##
        # if player['qmprog'] in ("g03", "g09", "g16"):

        ##
        ##  Molcas block
        ##
        # if player['qmprog'] == "molcas":

        # else:  ## Only relax the charge distribution

        #     if internal.player.qmprog in ("g03", "g09", "g16"):

        #         if cycle > 1:
        #             src = (
        #                 "simfiles"
        #                 + os.sep
        #                 + "step{:02d}".format(cycle - 1)
        #                 + os.sep
        #                 + "qm"
        #                 + os.sep
        #                 + "asec.chk"
        #             )
        #             dst = (
        #                 "simfiles"
        #                 + os.sep
        #                 + "step{:02d}".format(cycle)
        #                 + os.sep
        #                 + "qm"
        #                 + os.sep
        #                 + "asec.chk"
        #             )
        #             shutil.copyfile(src, dst)

        #         # internal.gaussian.make_charge_input(cycle, asec_charges)
        #         internal.gaussian.run_gaussian(cycle, "charge", internal.outfile)

        #         file = (
        #             "simfiles"
        #             + os.sep
        #             + "step{:02d}".format(cycle)
        #             + os.sep
        #             + "qm"
        #             + os.sep
        #             + "asec2.log"
        #         )
        #         internal.read_charges(file)

        #         ## Print new info for molecule[0]
        #         internal.outfile.write("\nNew values for molecule type 1:\n\n")
        #         internal.system.molecule[0].print_mol_info()

        # 			#if player['qmprog'] == "molcas"

        return position

    def reset(self):

        del self.step
