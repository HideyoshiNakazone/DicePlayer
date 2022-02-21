from diceplayer.DPpack.MolHandling import *
from diceplayer.DPpack.PTable import *
from diceplayer.DPpack.Misc import *

from diceplayer.DPpack.external.Dice import *

from typing import IO, Tuple, List, TextIO, Union

from numpy.core.numeric import partition
from numpy import random

import setproctitle
import subprocess
import os
import sys
import shutil
import textwrap
import types


dice_end_flag = "End of simulation"  # The normal end flag
dice_flag_line = -2  # must be in the line before the last
umaAng3_to_gcm3 = 1.6605  # Conversion between uma/Ang3 to g/cm3

max_seed = 4294967295  # Maximum allowed value for a seed (numpy)


class Internal:

    def __init__(self, infile: TextIO, outfile: TextIO) -> None:

        self.infile = infile
        self.outfile = outfile

        self.continued: bool = None

        self.system = System()

        self.player = self.Player()
        self.player_keywords = [a for a in dir(self.player) if not a.startswith(
            '__') and not callable(getattr(self.player, a))]

        self.dice = Dice(infile, outfile)
        self.dice_keywords = [a for a in dir(self.dice) if not a.startswith(
            '__') and not callable(getattr(self.dice, a))]

        self.gaussian = self.Gaussian()
        self.gaussian_keywords = [a for a in dir(self.gaussian) if not a.startswith(
            '__') and not callable(getattr(self.gaussian, a))]

        # self.molcas = self.Molcas()
        # self.molcas_keywords = [a for a in dir(self.molcas) if not a.startswith('__') and not callable(getattr(self.molcas, a))]

        # Constanst that shall be set for global use

        self.tol_rms_force = 3e-4		# Hartree/Bohr
        self.tol_max_force = 4.5e-4		# Hartree/Bohr
        self.tol_rms_step = 1.2e-3		# Bohr
        self.tol_max_step = 1.8e-3		# Bohr
        self.trust_radius = None

        # Dice:
        self.combrule = None

    def read_keywords(self) -> None:

        try:
            controlfile = self.infile.readlines()
        except EnvironmentError:
            sys.exit("Error: cannot read file {}".format(self.infile))

        for line in controlfile:

            key, value = line.partition("=")[::2]  # Discards the '='
            key = key.strip().lower()
            if key in ('title', 'keywords'):
                value = value.strip()
            else:
                value = value.split()

            # Read the Diceplayer related keywords
            # 'value' is not empty!
            if key in self.player_keywords and len(value) != 0:

                if key == 'qmprog' and value[0].lower() in ("g03", "g09", "g16", "molcas"):

                    setattr(self.player, key, value[0].lower())

                    if self.player.qmprog in ("g03", "g09", "g16"):

                        self.gaussian.qmprog = self.player.qmprog

                    # if self.player.qmprog == "molcas":

                    # 	pass

                elif key == 'opt' and value[0].lower() in ("yes", "no", "ts"):
                    setattr(self.player, key, value[0].lower())

                # elif key == 'zipprog' and value[0].lower() in ("zip", "gzip", "bzip"):
                    #player[key] = value[0].lower()

                elif key in ('lps', 'ghosts') and value[0].lower() in ("yes", "no"):
                    setattr(self.player, key, value[0].lower())

                elif key in ('readhessian', 'vdwforces') and value[0].lower() in ("yes", "no"):
                    setattr(self.player, key, value[0].lower())

                elif key in ('maxcyc', 'nprocs', 'altsteps', 'switchcyc'):
                    err = "Error: expected a positive integer for keyword {} in file {}".format(
                        key, self.infile)
                    try:
                        new_value = int(value[0])
                        if new_value >= 1:
                            setattr(self.player, key, new_value)
                        elif key == 'altsteps' and new_value == 0:
                            setattr(self.player, key, 0)
                    except ValueError:
                        sys.exit(err)

                elif key == 'maxstep':  # Cannot be less than 0.01
                    err = "Error: expected a float greater than 0.01 for keyword {} in file {}".format(
                        key, self.infile)
                    try:
                        new_value = float(value[0])
                        if new_value < 0.01:
                            sys.exit(err)
                        else:
                            setattr(self.player, key, new_value)
                    except ValueError:
                        sys.exit(err)

            # Read the Dice related keywords
            # 'value' is not empty!
            elif key in self.dice_keywords and len(value) != 0:

                if key == 'title':
                    setattr(self.dice, key, value)

                elif key in ('ljname', 'outname', 'progname'):
                    setattr(self.dice, key, value[0])

                elif key == 'randominit':
                    if value in ('always', 'first'):
                        setattr(self.dice, key, value[0])

                elif key in ('ncores', 'isave'):
                    err = "Error: expected a positive integer for keyword {} in file {}".format(
                        key, self.infile)
                    if not value[0].isdigit():
                        sys.exit(err)
                    new_value = int(value[0])
                    if new_value >= 1:
                        setattr(self.dice, key, new_value)

                elif key in ('temp', 'press', 'dens'):  # Cannot be less than 1e-10
                    err = "Error: expected a positive float for keyword {} in file {}".format(
                        key, self.infile)
                    try:
                        new_value = float(value[0])
                        if new_value < 1e-10:
                            sys.exit(err)
                        else:
                            setattr(self.dice, key, new_value)
                    except ValueError:
                        sys.exit(err)

                # If defined, must be well defined (only positive integer values)
                elif key == 'nmol':
                    err = "Error: expected 1 to 4 positive integers for keyword {} in file {}".format(
                        key, self.infile)
                    args = min(4, len(value))
                    for i in range(args):
                        if value[i].isdigit():
                            new_value = int(value[i])
                            if new_value < 1:
                                sys.exit(err)
                            else:
                                getattr(self.dice, key).append(new_value)
                        elif i == 0:
                            sys.exit(err)
                        else:
                            break

                # If defined, must be well defined (only positive integer values)
                elif key == 'nstep':
                    err = "Error: expected 2 or 3 positive integers for keyword {} in file {}".format(
                        key, self.infile)
                    if len(value) < 2:
                        sys.exit(err)
                    args = min(3, len(value))
                    for i in range(args):
                        if value[i].isdigit():
                            new_value = int(value[i])
                            if new_value < 1:
                                sys.exit(err)
                            else:
                                getattr(self.dice, key).append(new_value)
                        elif i < 2:
                            sys.exit(err)
                        else:
                            break

            # Read the Gaussian related keywords
            # 'value' is not empty!
            elif key in self.gaussian_keywords and len(value) != 0:

                if key == 'mem':  # Memory in MB (minimum of 100)
                    err = "Error: expected a positive integer for keyword {} in file {}".format(
                        key, self.infile)
                    if not value[0].isdigit():
                        sys.exit(err)
                    new_value = int(value[0])
                    if new_value >= 100:
                        setattr(self.gaussian, key, new_value)

                elif key == 'keywords':
                    setattr(self.gaussian, key, value)

                # If defined, must be well defined (2 integer values)
                elif key == 'chgmult':
                    err = "Error: expected 2 integers for keyword {} in file {}".format(
                        key, self.infile)
                    if len(value) < 2:
                        sys.exit(err)
                    for i in range(2):
                        try:
                            getattr(self.gaussian, key)[i] = int(value[i])
                        except ValueError:
                            sys.exit(err)

                elif key == 'level':
                    setattr(self.gaussian, key, value[0])

                elif key in ('gmiddle', 'gbottom'):
                    setattr(self.gaussian, key, value[0])

                elif key == 'pop' and value[0].lower() in ("chelpg", "mk", "nbo"):
                    setattr(self.gaussian, key, value[0].lower())

            # ####  Read the Molcas related keywords
            # elif key in self.molcas_keywords and len(value) != 0:  ##  'value' is not empty!

            # 	if key == 'root': # If defined, must be well defined (only positive integer values)
            # 		err = "Error: expected a positive integer for keyword {} in file {}".format(key, self.infile)
            # 		if not value[0].isdigit():
            # 			sys.exit(err)
            # 		new_value = int(value[0])
            # 		if new_value >= 1:
            # 			setattr(self.molcas, key, new_value)

            # 	elif key in ('mbottom', 'orbfile'):
            # 		setattr(self.molcas, key, value[0])

            # 	elif key == 'basis':
            # 		setattr(self.molcas ,key, value[0])

            # #### End

    def check_keywords(self) -> None:

        min_steps = 20000

        if self.dice.ljname == None:
            sys.exit(
                "Error: 'ljname' keyword not specified in file {}".format(self.infile))

        if self.dice.outname == None:
            sys.exit(
                "Error: 'outname' keyword not specified in file {}".format(self.infile))

        if self.dice.dens == None:
            sys.exit(
                "Error: 'dens' keyword not specified in file {}".format(self.infile))

        if self.dice.nmol == 0:
            sys.exit("Error: 'nmol' keyword not defined appropriately in file {}".format(
                self.infile))

        if self.dice.nstep == 0:
            sys.exit("Error: 'nstep' keyword not defined appropriately in file {}".format(
                self.infile))

        # Check only if QM program is Gaussian:
        if self.player.qmprog in ("g03", "g09", "g16"):

            if self.gaussian.level == None:
                sys.exit(
                    "Error: 'level' keyword not specified in file {}".format(self.infile))

            if self.gaussian.gmiddle != None:
                if not os.path.isfile(self.gaussian.gmiddle):
                    sys.exit("Error: file {} not found".format(
                        self.gaussian.gmiddle))

            if self.gaussian.gbottom != None:
                if not os.path.isfile(self.gaussian.gbottom):
                    sys.exit("Error: file {} not found".format(
                        self.gaussian.gbottom))

            if self.gaussian.pop != "chelpg" and (self.player.ghosts == "yes" or self.player.lps == "yes"):
                sys.exit(
                    "Error: ghost atoms or lone pairs only available with 'pop = chelpg')")

        # Check only if QM program is Molcas:
        # if self.player.qmprog == "molcas":

        # 	if self.molcas.mbottom == None:
        # 		sys.exit("Error: 'mbottom' keyword not specified in file {}".format(self.infile))
        # 	else:
        # 		if not os.path.isfile(self.molcas.mbottom):
        # 			sys.exit("Error: file {} not found".format(self.molcas.mbottom))

        # 	if self.molcas.basis == None:
        # 		sys.exit("Error: 'basis' keyword not specified in file {}".format(self.infile))

        if self.player.altsteps != 0:

            # Verifica se tem mais de 1 molecula QM
            # (No futuro usar o RMSD fit para poder substituir todas as moleculas QM
            # no arquivo outname.xy - Need to change the __make_init_file!!)
            if self.dice.nmol[0] > 1:
                sys.exit(
                    "Error: altsteps > 0 only possible with 1 QM molecule (nmol = 1 n2 n3 n4)")

            # if not zero, altsteps cannot be less than min_steps
            self.player.altsteps = max(min_steps, self.player.altsteps)
            # altsteps value is always the nearest multiple of 1000
            self.player.altsteps = round(self.player.altsteps / 1000) * 1000

        for i in range(len(self.dice.nstep)):
            # nstep can never be less than min_steps
            self.dice.nstep[i] = max(min_steps, self.dice.nstep[i])
            # nstep values are always the nearest multiple of 1000
            self.dice.nstep[i] = round(self.dice.nstep[i] / 1000) * 1000

        # isave must be between 100 and 2000
        self.dice.isave = max(100, self.dice.isave)
        self.dice.isave = min(2000, self.dice.isave)
        # isave value is always the nearest multiple of 100
        self.dice.isave = round(self.dice.isave / 100) * 100

    def print_keywords(self) -> None:

        self.outfile.write("##########################################################################################\n"
                           "#############               Welcome to DICEPLAYER version 1.0                #############\n"
                           "##########################################################################################\n"
                           "\n")
        self.outfile.write("Your python version is {}\n".format(sys.version))
        self.outfile.write("\n")
        self.outfile.write(
            "Program started on {}\n".format(weekday_date_time()))
        self.outfile.write("\n")
        self.outfile.write("Environment variables:\n")
        for var in env:
            self.outfile.write("{} = {}\n".format(var,
                                                  (os.environ[var] if var in os.environ else "Not set")))

        self.outfile.write("\n==========================================================================================\n"
                           "                         CONTROL variables being used in this run:\n"
                           "------------------------------------------------------------------------------------------\n"
                           "\n")

        for key in sorted(self.player_keywords):
            if getattr(self.player, key) != None:
                if isinstance(getattr(self.player, key), list):
                    string = " ".join(str(x)
                                      for x in getattr(self.player, key))
                    self.outfile.write("{} = {}\n".format(key, string))
                else:
                    self.outfile.write("{} = {}\n".format(
                        key, getattr(self.player, key)))

        self.outfile.write("\n")

        self.outfile.write("------------------------------------------------------------------------------------------\n"
                           "                         DICE variables being used in this run:\n"
                           "------------------------------------------------------------------------------------------\n"
                           "\n")

        for key in sorted(self.dice_keywords):
            if getattr(self.dice, key) != None:
                if isinstance(getattr(self.dice, key), list):
                    string = " ".join(str(x) for x in getattr(self.dice, key))
                    self.outfile.write("{} = {}\n".format(key, string))
                else:
                    self.outfile.write("{} = {}\n".format(
                        key, getattr(self.dice, key)))

        self.outfile.write("\n")

        if self.player.qmprog in ("g03", "g09", "g16"):

            self.outfile.write("------------------------------------------------------------------------------------------\n"
                               "                         GAUSSIAN variables being used in this run:\n"
                               "------------------------------------------------------------------------------------------\n"
                               "\n")

            for key in sorted(self.gaussian_keywords):
                if getattr(self.gaussian, key) != None:
                    if isinstance(getattr(self.gaussian, key), list):
                        string = " ".join(str(x)
                                          for x in getattr(self.gaussian, key))
                        self.outfile.write("{} = {}\n".format(key, string))
                    else:
                        self.outfile.write("{} = {}\n".format(
                            key, getattr(self.gaussian, key)))

            self.outfile.write("\n")

        # elif self.player.qmprog == "molcas":

        # 	self.outfile.write("------------------------------------------------------------------------------------------\n"
        # 			"                         MOLCAS variables being used in this run:\n"
        # 			"------------------------------------------------------------------------------------------\n"
        # 			"\n")

        # 	for key in sorted(molcas):
        # 		if molcas[key] != None:
        # 			if isinstance(molcas[key], list):
        # 				string = " ".join(str(x) for x in molcas[key])
        # 				self.outfile.write("{} = {}\n".format(key, string))
        # 			else:
        # 				self.outfile.write("{} = {}\n".format(key, molcas[key]))

        # 	self.outfile.write("\n")

    def read_potential(self) -> None:  # Deve ser atualizado para o uso de

        try:
            with open(self.dice.ljname) as file:
                ljfile = file.readlines()
        except EnvironmentError as err:
            sys.exit(err)

        combrule = ljfile.pop(0).split()[0]
        if combrule not in ("*", "+"):
            sys.exit(
                "Error: expected a '*' or a '+' sign in 1st line of file {}".format(self.dice.ljname))
        self.dice.combrule = combrule

        ntypes = ljfile.pop(0).split()[0]
        if not ntypes.isdigit():
            sys.exit("Error: expected an integer in the 2nd line of file {}".format(
                self.dice.ljname))
        ntypes = int(ntypes)

        if ntypes != len(self.dice.nmol):
            sys.exit("Error: number of molecule types in file {} must match that of 'nmol' keyword in file {}".format(
                self.dice.ljname, self.infile))
        line = 2
        for i in range(ntypes):

            line += 1
            nsites, molname = ljfile.pop(0).split()[:2]

            if not nsites.isdigit():
                sys.exit("Error: expected an integer in line {} of file {}".format(
                    line, self.dice.ljname))

            if molname is None:
                sys.exit("Error: expected a molecule name in line {} of file {}".format(
                    line, self.dice.ljname))

            nsites = int(nsites)

            self.system.add_type(nsites, Molecule(molname))

            for j in range(nsites):

                line += 1
                new_atom = ljfile.pop(0).split()

                if len(new_atom) < 8:
                    sys.exit("Error: expected at least 8 fields in line {} of file {}".format(
                        line, self.dice.ljname))

                if not new_atom[0].isdigit():
                    sys.exit("Error: expected an integer in field 1, line {} of file {}".format(
                        line, self.dice.ljname))
                lbl = int(new_atom[0])

                if not new_atom[1].isdigit():
                    sys.exit("Error: expected an integer in field 2, line {} of file {}".format(
                        line, self.dice.ljname))

                atnumber = int(new_atom[1])
                if atnumber == ghost_number and i == 0:  # Ghost atom not allowed in the QM molecule
                    sys.exit("Error: found a ghost atom in line {} of file {}".format(
                        line, self.dice.ljname))
                na = atnumber

                try:
                    rx = float(new_atom[2])
                except:
                    sys.exit("Error: expected a float in field 3, line {} of file {}".format(
                        line, self.dice.ljname))

                try:
                    ry = float(new_atom[3])
                except:
                    sys.exit("Error: expected a float in field 4, line {} of file {}".format(
                        line, self.dice.ljname))

                try:
                    rz = float(new_atom[4])
                except:
                    sys.exit("Error: expected a float in field 5, line {} of file {}".format(
                        line, self.dice.ljname))

                try:
                    chg = float(new_atom[5])
                except:
                    sys.exit("Error: expected a float in field 6, line {} of file {}".format(
                        line, self.dice.ljname))

                try:
                    eps = float(new_atom[6])
                except:
                    sys.exit("Error: expected a float in field 7, line {} of file {}".format(
                        line, self.dice.ljname))

                try:
                    sig = float(new_atom[7])
                except:
                    sys.exit("Error: expected a float in field 8, line {} of file {}".format(
                        line, self.dice.ljname))

                mass = atommass[na]

                if len(new_atom) > 8:
                    masskey, mass = new_atom[8].partition("=")[::2]
                    if masskey.lower() == 'mass' and len(mass) != 0:
                        try:
                            new_mass = float(mass)
                            if new_mass > 0:
                                mass = new_mass
                        except:
                            sys.exit(
                                "Error: expected a positive float after 'mass=' in field 9, line {} of file {}".format(
                                    line, self.dice.ljname))

                self.system.molecule[i].add_atom(
                    Atom(lbl, na, rx, ry, rz, chg, eps, sig))

        to_delete = ['lbl', 'na', 'rx', 'ry',
                     'rz', 'chg', 'eps', 'sig', 'mass']
        for _var in to_delete:
            if _var in locals() or _var in globals():
                exec(f'del {_var}')

    def print_potential(self) -> None:

        formatstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f} {:>9.4f}\n"
        self.outfile.write("\n"
                           "==========================================================================================\n")
        self.outfile.write(
            "                    Potential parameters from file {}:\n".format(self.dice.ljname))
        self.outfile.write("------------------------------------------------------------------------------------------\n"
                           "\n")

        self.outfile.write("Combination rule: {}\n".format(self.dice.combrule))
        self.outfile.write("Types of molecules: {}\n\n".format(
            len(self.system.molecule)))

        i = 0
        for mol in self.system.molecule:
            i += 1
            self.outfile.write(
                "{} atoms in molecule type {}:\n".format(len(mol.atom), i))
            self.outfile.write("---------------------------------------------------------------------------------\n"
                               "Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass\n")
            self.outfile.write(
                "---------------------------------------------------------------------------------\n")

            for atom in mol.atom:

                self.outfile.write(formatstr.format(atom.lbl, atom.na, atom.rx, atom.ry, atom.rz,
                                                    atom.chg, atom.eps, atom.sig, atom.mass))

            self.outfile.write("\n")

        if self.player.ghosts == "yes" or self.player.lps == "yes":
            self.outfile.write("\n"
                               "------------------------------------------------------------------------------------------\n"
                               "                    Aditional potential parameters:\n"
                               "------------------------------------------------------------------------------------------\n")

        # if player['ghosts'] == "yes":

        # 	self.outfile.write("\n")
        # 	self.outfile.write("{} ghost atoms appended to molecule type 1 at:\n".format(len(ghost_types)))
        # 	self.outfile.write("---------------------------------------------------------------------------------\n")

        # 	atoms_string = ""
        # 	for ghost in ghost_types:
        # 		for atom in ghost['numbers']:
        # 			atom_sym = atomsymb[ molecules[0][atom - 1]['na'] ].strip()
        # 			atoms_string += "{}{} ".format(atom_sym,atom)

        # 		if ghost['type'] == "g":
        # 			self.outfile.write(textwrap.fill("* Geometric center of atoms {}".format(atoms_string), 80))
        # 		elif ghost['type'] == "m":
        # 			self.outfile.write(textwrap.fill("* Center of mass of atoms {}".format(atoms_string), 80))
        # 		elif ghost['type'] == "z":
        # 			self.outfile.write(textwrap.fill("* Center of atomic number of atoms {}".format(atoms_string), 80))

        # 		self.outfile.write("\n")

        # if player['lps'] == 'yes':

        # 	self.outfile.write("\n")
        # 	self.outfile.write("{} lone pairs appended to molecule type 1:\n".format(len(lp_types)))
        # 	self.outfile.write("---------------------------------------------------------------------------------\n")

        # 	for lp in lp_types:
        # 		# LP type 1 or 2
        # 		if lp['type'] in (1, 2):
        # 			atom1_num = lp['numbers'][0]
        # 			atom1_sym = atomsymb[ molecules[0][atom1_num - 1]['na'] ].strip()
        # 			atom2_num = lp['numbers'][1]
        # 			atom2_sym = atomsymb[ molecules[0][atom2_num - 1]['na'] ].strip()
        # 			atom3_num = lp['numbers'][2]
        # 			atom3_sym = atomsymb[ molecules[0][atom3_num - 1]['na'] ].strip()

        # 			self.outfile.write(textwrap.fill(
        # 			"* Type {} on atom {}{} with {}{} {}{}. Alpha = {:<5.1f} Deg and D = {:<4.2f} Angs".format(
        # 				lp['type'], atom1_sym, atom1_num, atom2_sym, atom2_num, atom3_sym, atom3_num, lp['alpha'],
        # 				lp['dist']), 86))
        # 			self.outfile.write("\n")

        # 		# Other LP types

        self.outfile.write("\n"
                           "==========================================================================================\n")

    def check_executables(self) -> None:

        self.outfile.write("\n")
        self.outfile.write(90 * "=")
        self.outfile.write("\n\n")

        dice_path = shutil.which(self.dice.progname)
        if dice_path != None:
            self.outfile.write("Program {} found at {}\n".format(
                self.dice.progname, dice_path))
            self.dice.path = dice_path
        else:
            sys.exit("Error: cannot find dice executable")

        qmprog_path = shutil.which(self.gaussian.qmprog)
        if qmprog_path != None:
            self.outfile.write("Program {} found at {}\n".format(
                self.gaussian.qmprog, qmprog_path))
            self.gaussian.path = qmprog_path
        else:
            sys.exit("Error: cannot find {} executable".format(
                self.gaussian.qmprog))

        if self.gaussian.qmprog in ("g03", "g09", "g16"):
            formchk_path = shutil.which("formchk")
            if formchk_path != None:
                self.outfile.write(
                    "Program formchk found at {}\n".format(formchk_path))
            else:
                sys.exit("Error: cannot find formchk executable")

    def dice_start(self, cycle: int):

        self.dice.configure(self.player.initcyc,
                            self.player.nprocs,
                            self.player.altsteps,
                            self.system.nmols,
                            self.system.molecule)

        self.dice.start(cycle)

        self.dice.reset()

    def calculate_step(self, cycle: int, gradient: np.ndarray,
                       hessian: np.ndarray) -> np.ndarray:

        invhessian = np.linalg.inv(hessian)
        pre_step = -1 * np.matmul(invhessian, gradient.T).T
        maxstep = np.amax(np.absolute(pre_step))
        factor = min(1, self.player.maxstep/maxstep)
        step = factor * pre_step

        self.outfile.write("\nCalculated step-{}:\n".format(cycle))
        pre_step_list = pre_step.tolist()

        self.outfile.write("-----------------------------------------------------------------------\n"
                           "Center     Atomic                          Step (Bohr)\n"
                           "Number     Number              X                Y                Z\n"
                           "-----------------------------------------------------------------------\n")
        for i in range(len(self.system.molecule[0].atom)):
            self.outfile.write("  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
                i + 1, self.system.molecule[0].atom[i].na,
                pre_step_list.pop(0), pre_step_list.pop(0), pre_step_list.pop(0)))

        self.outfile.write(
            "-----------------------------------------------------------------------\n")

        self.outfile.write("Maximum step is {:>11.6}\n".format(maxstep))
        self.outfile.write("Scaling factor = {:>6.4f}\n".format(factor))
        self.outfile.write("\nFinal step (Bohr):\n")
        step_list = step.tolist()

        self.outfile.write("-----------------------------------------------------------------------\n"
                           "Center     Atomic                          Step (Bohr)\n"
                           "Number     Number              X                Y                Z\n"
                           "-----------------------------------------------------------------------\n")
        for i in range(len(self.system.molecule[0].atom)):
            self.outfile.write("  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
                i + 1, self.system.molecule[0].atom[i].na,
                step_list.pop(0), step_list.pop(0), step_list.pop(0)))

        self.outfile.write(
            "-----------------------------------------------------------------------\n")

        step_max = np.amax(np.absolute(step))
        step_rms = np.sqrt(np.mean(np.square(step)))

        self.outfile.write("  Max Step = {:>14.9f}      RMS Step = {:>14.9f}\n\n".format(
            step_max, step_rms))

        return step

    # I still have to talk with Herbet about this function
    def populate_asec_vdw(self, cycle):

        # Both asec_charges and vdw_meanfield will utilize the Molecule() class and Atoms() with some None elements

        asec_charges = Molecule()  	# (lbl=None, na=None, rx, ry, rz, chg, eps=None, sig=None)
        vdw_meanfield = Molecule() 	# (lbl=None, na=None, rx, ry, rz, chg=None, eps, sig)

        if self.dice.nstep[-1] % self.dice.isave == 0:
            nconfigs = round(self.dice.nstep[-1] / self.dice.isave)
        else:
            nconfigs = int(self.dice.nstep[-1] / self.dice.isave)

        norm_factor = nconfigs * self.player.nprocs

        nsitesref = len(self.system.molecule[0].atom) + len(
            self.system.molecule[0].ghost_atoms) + len(self.system.molecule[0].lp_atoms)

        nsites_total = self.dice.nmol[0] * nsitesref
        for i in range(1, len(self.dice.nmol)):
            nsites_total += self.dice.nmol[i] * \
                len(self.system.molecule[i].atom)

        thickness = []
        picked_mols = []

        for proc in range(1, self.player.nprocs + 1):  # Run over folders

            simdir = "simfiles"
            path = simdir + os.sep + \
                "step{:02d}".format(cycle) + os.sep + "p{:02d}".format(proc)
            file = path + os.sep + self.dice.outname + ".xyz"
            if not os.path.isfile(file):
                sys.exit("Error: cannot find file {}".format(file))
            try:
                with open(file) as xyzfh:
                    xyzfile = xyzfh.readlines()
            except:
                sys.exit("Error: cannot open file {}".format(file))

            for config in range(nconfigs):  # Run over configs in a folder

                if int(xyzfile.pop(0).split()[0]) != nsites_total:
                    sys.exit(
                        "Error: wrong number of sites in file {}".format(file))

                box = xyzfile.pop(0).split()[-3:]
                box = [float(box[0]), float(box[1]), float(box[2])]
                sizes = self.system.molecule[0].sizes_of_molecule()
                thickness.append(min([(box[0] - sizes[0])/2, (box[1] - sizes[1])/2,
                                      (box[2] - sizes[2])/2]))

                # Skip the first (reference) molecule
                xyzfile = xyzfile[nsitesref:]
                mol_count = 0
                for type in range(len(self.dice.nmol)):  # Run over types of molecules

                    if type == 0:
                        nmols = self.dice.nmol[0] - 1
                    else:
                        nmols = self.dice.nmol[type]

                    for mol in range(nmols):  # Run over molecules of each type

                        new_molecule = Molecule(
                            self.system.molecule[type].molnale)
                        # Run over sites of each molecule
                        for site in range(len(self.system.molecule[types].atom)):

                            new_molecule.append({})
                            line = xyzfile.pop(0).split()

                            if line[0].title() != atomsymb[self.system.molecule[type].atom[site].na.strip()]:
                                sys.exit("Error reading file {}".format(file))

                            new_molecule.add_atom(Atom(self.system.molecule[type].atom[site].lbl,
                                                       self.system.molecule[type].atom[site].na,
                                                       self.system.molecule[type].atom[site].float(
                                                           line[1]),
                                                       self.system.molecule[type].atom[site].float(
                                                           line[2]),
                                                       self.system.molecule[type].atom[site].float(
                                                           line[3]),
                                                       self.system.molecule[type].atom[site].chg,
                                                       self.system.molecule[type].atom[site].eps,
                                                       self.system.molecule[type].atom[site].sig))

                        dist = self.system.molecule[0].minimum_distance(
                            new_molecule)
                        if dist < thickness[-1]:
                            mol_count += 1
                            for atom in new_molecule:
                                asec_charges.append({})
                                vdw_meanfield.append({})

                                asec_charges[-1]['rx'] = atom['rx']
                                asec_charges[-1]['ry'] = atom['ry']
                                asec_charges[-1]['rz'] = atom['rz']
                                asec_charges[-1]['chg'] = atom['chg'] / \
                                    norm_factor

                                if self.player.vdwforces == "yes":
                                    vdw_meanfield[-1]['rx'] = atom['rx']
                                    vdw_meanfield[-1]['ry'] = atom['ry']
                                    vdw_meanfield[-1]['rz'] = atom['rz']
                                    vdw_meanfield[-1]['eps'] = atom['eps']
                                    vdw_meanfield[-1]['sig'] = atom['sig']

                        # ####  Read lines with ghosts or lps in molecules of type 0 (reference)
                        # ####  and, if dist < thickness, appends to asec
                        # if type == 0:
                        # 	for ghost in ghost_atoms:
                        # 		line = xyzfile.pop(0).split()
                        # 		if line[0] != dice_ghost_label:
                        # 			sys.exit("Error reading file {}".format(file))
                        # 		if dist < thickness[-1]:
                        # 			asec_charges.append({})
                        # 			asec_charges[-1]['rx'] = float(line[1])
                        # 			asec_charges[-1]['ry'] = float(line[2])
                        # 			asec_charges[-1]['rz'] = float(line[3])
                        # 			asec_charges[-1]['chg'] = ghost['chg'] / norm_factor

                        # 	for lp in lp_atoms:
                        # 		line = xyzfile.pop(0).split()
                        # 		if line[0] != dice_ghost_label:
                        # 			sys.exit("Error reading file {}".format(file))
                        # 		if dist < thickness[-1]:
                        # 			asec_charges.append({})
                        # 			asec_charges[-1]['rx'] = float(line[1])
                        # 			asec_charges[-1]['ry'] = float(line[2])
                        # 			asec_charges[-1]['rz'] = float(line[3])
                        # 			asec_charges[-1]['chg'] = lp['chg'] / norm_factor

                picked_mols.append(mol_count)

        self.outfile.write("Done\n")

        string = "In average, {:^7.2f} molecules ".format(
            sum(picked_mols)/norm_factor)
        string += "were selected from each of the {} configurations ".format(
            len(picked_mols))
        string += "of the production simulations to form the ASEC, comprising a shell with "
        string += "minimum thickness of {:>6.2f} Angstrom\n".format(
            sum(thickness)/norm_factor)

        self.outfile.write(textwrap.fill(string, 86))
        self.outfile.write("\n")

        otherfh = open("ASEC.dat", "w")
        for charge in asec_charges:
            otherfh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
                charge['rx'], charge['ry'], charge['rz'], charge['chg']))
        otherfh.close()

        return asec_charges

    # ----------------------------------------------------------------------------------------------------------------------------------------
    # Gaussian related methods
    # ----------------------------------------------------------------------------------------------------------------------------------------
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

        degrees = 3 * len(self.system.molecule[0].atom)
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
        fh.write("-----------------------------------------------------------------------\n"
                 "Center     Atomic                     Forces (Hartree/Bohr)\n"
                 "Number     Number              X                Y                Z\n"
                 "-----------------------------------------------------------------------\n")
        for i in range(len(self.system.molecule[0].atom)):
            fh.write("  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
                i + 1, self.system.molecule[0].atom[i].na, forces.pop(0), forces.pop(0), forces.pop(0)))

        fh.write(
            "-----------------------------------------------------------------------\n")

        force_max = np.amax(np.absolute(gradient))
        force_rms = np.sqrt(np.mean(np.square(gradient)))

        fh.write("  Max Force = {:>14.9f}      RMS Force = {:>14.9f}\n\n".format(
            force_max, force_rms))

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

        degrees = 3 * len(self.system.molecule[0].atom)
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

        degrees = 3 * len(self.system.molecule[0].atom)
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

    def print_grad_hessian(self, cycle: int,
                           cur_gradient: np.ndarray, hessian: np.ndarray
                           ) -> None:

        try:
            fh = open("grad_hessian.dat", "w")
        except:
            sys.exit("Error: cannot open file grad_hessian.dat")

        fh.write("Optimization cycle: {}\n".format(cycle))
        fh.write("Cartesian Gradient\n")
        degrees = 3 * len(self.system.molecule[0].atom)
        for i in range(degrees):
            fh.write(" {:>11.8g}".format(cur_gradient[i]))
            if (i + 1) % 5 == 0 or i == degrees - 1:
                fh.write("\n")

        fh.write("Cartesian Force Constants\n")
        n = int(np.sqrt(2*degrees))
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
        if self.gaussian.mem != None:
            fh.write("%Mem={}MB\n".format(self.gaussian.mem))
        fh.write("%Nprocs={}\n".format(self.player.nprocs * self.dice.ncores))

        kword_line = "#P " + str(self.gaussian.level)

        if self.gaussian.keywords != None:
            kword_line += " " + self.gaussian.keywords

        if self.player.opt == 'yes':
            kword_line += " Force"

        # kword_line += " Charge"
        kword_line += " NoSymm"
        kword_line += " Pop={} Density=Current".format(self.gaussian.pop)

        if cycle > 1:
            kword_line += " Guess=Read"

        fh.write(textwrap.fill(kword_line, 90))
        fh.write("\n")

        fh.write("\nForce calculation - Cycle number {}\n".format(cycle))
        fh.write("\n")
        fh.write("{},{}\n".format(
            self.gaussian.chgmult[0], self.gaussian.chgmult[1]))

        for atom in self.system.molecule[0].atom:
            symbol = atomsymb[atom.na]
            fh.write("{:<2s}    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(symbol,
                                                                            atom.rx, atom.ry, atom.rz))

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
        # if self.gaussian.gmiddle != None:
        # 	if not os.path.isfile(self.gaussian.gmiddle):
        # 		sys.exit("Error: cannot find file {} in main directory".format(
        # 																self.gaussian.gmiddle))
        # 	try:
        # 		with open(self.gaussian.gmiddle) as gmiddlefile:
        # 			gmiddle = gmiddlefile.readlines()
        # 	except:
        # 		sys.exit("Error: cannot open file {}".format(self.gaussian.gmiddle))

        # 	for line in gmiddle:
        # 		fh.write(line)

        # 	fh.write("\n")

        # ## Write the ASEC:
        # for charge in asec_charges:
        # 	fh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
        # 							charge['rx'], charge['ry'], charge['rz'], charge['chg']))

        fh.write("\n")

        # ## If gbottom file was informed, write its contents in asec.gjf
        # if self.gaussian.gbottom != None:
        # 	if not os.path.isfile(self.gaussian.gbottom):
        # 		sys.exit("Error: cannot find file {} in main directory".format(
        # 																self.gaussian.gbottom))
        # 	try:
        # 		with open(self.gaussian.gbottom) as gbottomfile:
        # 			gbottom = gbottomfile.readlines()
        # 	except:
        # 		sys.exit("Error: cannot open file {}".format(self.gaussian.gbottom))

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
        for atom in self.system.molecule[0].atom:
            line = glogfile.pop(0).split()
            atom_str = line[1]
            charge = float(line[2])
            atom.chg = charge
            fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))

        # if self.gaussian.pop == "chelpg":
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

    class Player:

        def __init__(self) -> None:

            self.maxcyc = None
            self.nprocs = 1
            self.switchcyc = 3
            self.altsteps = 20000
            self.maxstep = .3
            self.opt = "yes"
            self.freq = "no"
            self.readhessian = "no"
            self.lps = "no"
            self.ghosts = "no"
            self.vdwforces = "no"
            self.tol_factor = 1.2
            self.qmprog = "g16"

            self.initcyc = 1

    class Gaussian:

        def __init__(self) -> None:

            self.qmprog = "g09"
            self.path = None

            self.mem = None
            self.keywords = None
            self.chgmult = [0, 1]
            self.gmiddle = None   # In each case, if a filename is given, its content will be
            self.gbottom = None   # inserted in the gaussian input
            self.pop = "chelpg"
            self.level = None

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

            fh.write("\nCalculation of {}s initiated with Gaussian on {}\n".format(
                type, date_time()))

            if shutil.which("bash") != None:
                exit_status = subprocess.call(
                    ["bash", "-c", "exec -a {}-step{} {} {}".format(self.qmprog, cycle, self.qmprog, infile)])
            else:
                exit_status = subprocess.call([self.qmprog, infile])

            if exit_status != 0:
                sys.exit("Gaussian process did not exit properly")

            fh.write("Calculation of {}s finished on {}\n".format(
                type, date_time()))

            os.chdir(work_dir)

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

    # class Molcas:

    # 	def __init(self):

    # 		self.orbfile = "input.exporb"
    # 		self.root = 1

    # 		self.mbottom = None
    # 		self.basis = None
