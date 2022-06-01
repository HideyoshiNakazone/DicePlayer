import os
import shutil
import sys
import textwrap
import types
from typing import TextIO

import yaml

from diceplayer.DPpack.Environment.Atom import Atom
from diceplayer.DPpack.Environment.Molecule import Molecule
from diceplayer.DPpack.Environment.System import System
from diceplayer.DPpack.External.Dice import Dice
from diceplayer.DPpack.External.Gaussian import Gaussian
from diceplayer.DPpack.Utils.Misc import *
from diceplayer.DPpack.Utils.PTable import *
from diceplayer.DPpack.Utils.StepDTO import StepDTO
from diceplayer.DPpack.Utils.Validations import NotNull

env = ["OMP_STACKSIZE"]


class Player:

    maxcyc = None
    opt = None
    nprocs = None
    qmprog = None
    lps = None
    ghosts = None
    altsteps = None

    combrule = None

    TOL_RMS_FORCE = 3e-4
    TOL_MAX_FORCE = 4.5e-4
    TOL_RMS_STEP = 1.2e-3
    TOL_MAX_SET = 1.8e-3
    TRUST_RADIUS = None

    continued: bool = False

    def __init__(self, infile: TextIO, outfile: TextIO) -> None:

        self.infile = infile
        self.outfile = outfile

        self.system = System()

        self.dice = Dice(infile, outfile)
        self.dice_keywords = [
            a
            for a in dir(self.dice)
            if not a.startswith("__") and not callable(getattr(self.dice, a))
        ]

        self.gaussian = Gaussian()
        self.gaussian_keywords = [
            a
            for a in dir(self.gaussian)
            if not a.startswith("__") and not callable(getattr(self.gaussian, a))
        ]

    @NotNull(
        requiredArgs=["maxcyc", "opt", "nprocs", "qmprog", "lps", "ghosts", "altsteps"]
    )
    def updateKeywords(self, **data):
        self.__dict__.update(**data)

    def read_keywords(self) -> None:

        with self.infile as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)

        self.updateKeywords(**data.get("diceplayer"))
        self.dice.updateKeywords(**data.get("dice"))
        self.gaussian.updateKeywords(**data.get("gaussian"))

    def check_keywords(self) -> None:

        min_steps = 20000

        if self.dice.ljname == None:
            sys.exit(
                "Error: 'ljname' keyword not specified in file {}".format(self.infile)
            )

        if self.dice.outname == None:
            sys.exit(
                "Error: 'outname' keyword not specified in file {}".format(self.infile)
            )

        if self.dice.dens == None:
            sys.exit(
                "Error: 'dens' keyword not specified in file {}".format(self.infile)
            )

        if self.dice.nmol == 0:
            sys.exit(
                "Error: 'nmol' keyword not defined appropriately in file {}".format(
                    self.infile
                )
            )

        if self.dice.nstep == 0:
            sys.exit(
                "Error: 'nstep' keyword not defined appropriately in file {}".format(
                    self.infile
                )
            )

        # Check only if QM program is Gaussian:
        if self.qmprog in ("g03", "g09", "g16"):

            if self.gaussian.level == None:
                sys.exit(
                    "Error: 'level' keyword not specified in file {}".format(
                        self.infile
                    )
                )

            if self.gaussian.gmiddle != None:
                if not os.path.isfile(self.gaussian.gmiddle):
                    sys.exit("Error: file {} not found".format(self.gaussian.gmiddle))

            if self.gaussian.gbottom != None:
                if not os.path.isfile(self.gaussian.gbottom):
                    sys.exit("Error: file {} not found".format(self.gaussian.gbottom))

            if self.gaussian.pop != "chelpg" and (
                self.ghosts == "yes" or self.lps == "yes"
            ):
                sys.exit(
                    "Error: ghost atoms or lone pairs only available with 'pop = chelpg')"
                )

        # Check only if QM program is Molcas:
        # if self.qmprog == "molcas":

        # 	if self.molcas.mbottom == None:
        # 		sys.exit("Error: 'mbottom' keyword not specified in file {}".format(self.infile))
        # 	else:
        # 		if not os.path.isfile(self.molcas.mbottom):
        # 			sys.exit("Error: file {} not found".format(self.molcas.mbottom))

        # 	if self.molcas.basis == None:
        # 		sys.exit("Error: 'basis' keyword not specified in file {}".format(self.infile))

        if self.altsteps != 0:

            # Verifica se tem mais de 1 molecula QM
            # (No futuro usar o RMSD fit para poder substituir todas as moleculas QM
            # no arquivo outname.xy - Need to change the __make_init_file!!)
            if self.dice.nmol[0] > 1:
                sys.exit(
                    "Error: altsteps > 0 only possible with 1 QM molecule (nmol = 1 n2 n3 n4)"
                )

            # if not zero, altsteps cannot be less than min_steps
            self.altsteps = max(min_steps, self.altsteps)
            # altsteps value is always the nearest multiple of 1000
            self.altsteps = round(self.altsteps / 1000) * 1000

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

        self.outfile.write(
            "##########################################################################################\n"
            "#############               Welcome to DICEPLAYER version 1.0                #############\n"
            "##########################################################################################\n"
            "\n"
        )
        self.outfile.write("Your python version is {}\n".format(sys.version))
        self.outfile.write("\n")
        self.outfile.write("Program started on {}\n".format(weekday_date_time()))
        self.outfile.write("\n")
        self.outfile.write("Environment variables:\n")
        for var in env:
            self.outfile.write(
                "{} = {}\n".format(
                    var, (os.environ[var] if var in os.environ else "Not set")
                )
            )

        self.outfile.write(
            "\n==========================================================================================\n"
            "                         CONTROL variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
            "\n"
        )

        self.outfile.write("\n")

        self.outfile.write(
            "------------------------------------------------------------------------------------------\n"
            "                         DICE variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
            "\n"
        )

        for key in sorted(self.dice_keywords):
            if getattr(self.dice, key) != None:
                if isinstance(getattr(self.dice, key), list):
                    string = " ".join(str(x) for x in getattr(self.dice, key))
                    self.outfile.write("{} = {}\n".format(key, string))
                else:
                    self.outfile.write("{} = {}\n".format(key, getattr(self.dice, key)))

        self.outfile.write("\n")

        if self.qmprog in ("g03", "g09", "g16"):

            self.outfile.write(
                "------------------------------------------------------------------------------------------\n"
                "                         GAUSSIAN variables being used in this run:\n"
                "------------------------------------------------------------------------------------------\n"
                "\n"
            )

            for key in sorted(self.gaussian_keywords):
                if getattr(self.gaussian, key) != None:
                    if isinstance(getattr(self.gaussian, key), list):
                        string = " ".join(str(x) for x in getattr(self.gaussian, key))
                        self.outfile.write("{} = {}\n".format(key, string))
                    else:
                        self.outfile.write(
                            "{} = {}\n".format(key, getattr(self.gaussian, key))
                        )

            self.outfile.write("\n")

        # elif self.qmprog == "molcas":

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
                "Error: expected a '*' or a '+' sign in 1st line of file {}".format(
                    self.dice.ljname
                )
            )
        self.dice.combrule = combrule

        ntypes = ljfile.pop(0).split()[0]
        if not ntypes.isdigit():
            sys.exit(
                "Error: expected an integer in the 2nd line of file {}".format(
                    self.dice.ljname
                )
            )
        ntypes = int(ntypes)

        if ntypes != len(self.dice.nmol):
            sys.exit(
                "Error: number of molecule types in file {} must match that of 'nmol' keyword in file {}".format(
                    self.dice.ljname, self.infile
                )
            )
        line = 2
        for i in range(ntypes):

            line += 1
            nsites, molname = ljfile.pop(0).split()[:2]

            if not nsites.isdigit():
                sys.exit(
                    "Error: expected an integer in line {} of file {}".format(
                        line, self.dice.ljname
                    )
                )

            if molname is None:
                sys.exit(
                    "Error: expected a molecule name in line {} of file {}".format(
                        line, self.dice.ljname
                    )
                )

            nsites = int(nsites)

            self.system.add_type(nsites, Molecule(molname))

            for j in range(nsites):

                line += 1
                new_atom = ljfile.pop(0).split()

                if len(new_atom) < 8:
                    sys.exit(
                        "Error: expected at least 8 fields in line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                if not new_atom[0].isdigit():
                    sys.exit(
                        "Error: expected an integer in field 1, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )
                lbl = int(new_atom[0])

                if not new_atom[1].isdigit():
                    sys.exit(
                        "Error: expected an integer in field 2, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                atnumber = int(new_atom[1])
                if (
                    atnumber == ghost_number and i == 0
                ):  # Ghost atom not allowed in the QM molecule
                    sys.exit(
                        "Error: found a ghost atom in line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )
                na = atnumber

                try:
                    rx = float(new_atom[2])
                except:
                    sys.exit(
                        "Error: expected a float in field 3, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                try:
                    ry = float(new_atom[3])
                except:
                    sys.exit(
                        "Error: expected a float in field 4, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                try:
                    rz = float(new_atom[4])
                except:
                    sys.exit(
                        "Error: expected a float in field 5, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                try:
                    chg = float(new_atom[5])
                except:
                    sys.exit(
                        "Error: expected a float in field 6, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                try:
                    eps = float(new_atom[6])
                except:
                    sys.exit(
                        "Error: expected a float in field 7, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                try:
                    sig = float(new_atom[7])
                except:
                    sys.exit(
                        "Error: expected a float in field 8, line {} of file {}".format(
                            line, self.dice.ljname
                        )
                    )

                mass = atommass[na]

                if len(new_atom) > 8:
                    masskey, mass = new_atom[8].partition("=")[::2]
                    if masskey.lower() == "mass" and len(mass) != 0:
                        try:
                            new_mass = float(mass)
                            if new_mass > 0:
                                mass = new_mass
                        except:
                            sys.exit(
                                "Error: expected a positive float after 'mass=' in field 9, line {} of file {}".format(
                                    line, self.dice.ljname
                                )
                            )

                self.system.molecule[i].add_atom(
                    Atom(lbl, na, rx, ry, rz, chg, eps, sig)
                )

        to_delete = ["lbl", "na", "rx", "ry", "rz", "chg", "eps", "sig", "mass"]
        for _var in to_delete:
            if _var in locals() or _var in globals():
                exec(f"del {_var}")

    def print_potential(self) -> None:

        formatstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f} {:>9.4f}\n"
        self.outfile.write(
            "\n"
            "==========================================================================================\n"
        )
        self.outfile.write(
            "                    Potential parameters from file {}:\n".format(
                self.dice.ljname
            )
        )
        self.outfile.write(
            "------------------------------------------------------------------------------------------\n"
            "\n"
        )

        self.outfile.write("Combination rule: {}\n".format(self.dice.combrule))
        self.outfile.write(
            "Types of molecules: {}\n\n".format(len(self.system.molecule))
        )

        i = 0
        for mol in self.system.molecule:
            i += 1
            self.outfile.write(
                "{} atoms in molecule type {}:\n".format(len(mol.atom), i)
            )
            self.outfile.write(
                "---------------------------------------------------------------------------------\n"
                "Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass\n"
            )
            self.outfile.write(
                "---------------------------------------------------------------------------------\n"
            )

            for atom in mol.atom:

                self.outfile.write(
                    formatstr.format(
                        atom.lbl,
                        atom.na,
                        atom.rx,
                        atom.ry,
                        atom.rz,
                        atom.chg,
                        atom.eps,
                        atom.sig,
                        atom.mass,
                    )
                )

            self.outfile.write("\n")

        if self.ghosts == "yes" or self.lps == "yes":
            self.outfile.write(
                "\n"
                "------------------------------------------------------------------------------------------\n"
                "                    Aditional potential parameters:\n"
                "------------------------------------------------------------------------------------------\n"
            )

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

        self.outfile.write(
            "\n"
            "==========================================================================================\n"
        )

    def check_executables(self) -> None:

        self.outfile.write("\n")
        self.outfile.write(90 * "=")
        self.outfile.write("\n\n")

        dice_path = shutil.which(self.dice.progname)
        if dice_path != None:
            self.outfile.write(
                "Program {} found at {}\n".format(self.dice.progname, dice_path)
            )
            self.dice.path = dice_path
        else:
            sys.exit("Error: cannot find dice executable")

        qmprog_path = shutil.which(self.gaussian.qmprog)
        if qmprog_path != None:
            self.outfile.write(
                "Program {} found at {}\n".format(self.gaussian.qmprog, qmprog_path)
            )
            self.gaussian.path = qmprog_path
        else:
            sys.exit("Error: cannot find {} executable".format(self.gaussian.qmprog))

        if self.gaussian.qmprog in ("g03", "g09", "g16"):
            formchk_path = shutil.which("formchk")
            if formchk_path != None:
                self.outfile.write("Program formchk found at {}\n".format(formchk_path))
            else:
                sys.exit("Error: cannot find formchk executable")

    def dice_start(self, cycle: int):

        self.dice.configure(
            StepDTO(
                initcyc=self.initcyc,
                nprocs=self.nprocs,
                altsteps=self.altsteps,
                nmol=self.system.nmols,
                molecule=self.system.molecule,
            )
        )

        self.dice.start(cycle)

        self.dice.reset()

    def gaussian_start(self, cycle: int, geomsfh: TextIO):

        self.gaussian.configure(
            self.initcyc,
            self.nprocs,
            self.dice.ncores,
            self.altsteps,
            self.switchcyc,
            self.opt,
            self.system.nmols,
            self.system.molecule,
        )

        position = self.gaussian.start(cycle, self.outfile, self.readhessian)

        ## Update the geometry of the reference molecule
        self.system.update_molecule(position, self.outfile)

        ## Print new geometry in geoms.xyz
        self.system.print_geom(cycle, geomsfh)

        self.gaussian.reset()

    # I still have to talk with Herbet about this function
    def populate_asec_vdw(self, cycle):

        # Both asec_charges and vdw_meanfield will utilize the Molecule() class and Atoms() with some None elements

        asec_charges = Molecule(
            "ASEC_CHARGES"
        )  # (lbl=None, na=None, rx, ry, rz, chg, eps=None, sig=None)
        # vdw_meanfield = (
        #     Molecule()
        # )  # (lbl=None, na=None, rx, ry, rz, chg=None, eps, sig)

        if self.dice.nstep[-1] % self.dice.isave == 0:
            nconfigs = round(self.dice.nstep[-1] / self.dice.isave)
        else:
            nconfigs = int(self.dice.nstep[-1] / self.dice.isave)

        norm_factor = nconfigs * self.nprocs

        nsitesref = len(self.system.molecule[0].atom)
        # nsitesref = (
        #     len(self.system.molecule[0].atom)
        #     + len(self.system.molecule[0].ghost_atoms)
        #     + len(self.system.molecule[0].lp_atoms)
        # )

        nsites_total = self.dice.nmol[0] * nsitesref
        for i in range(1, len(self.dice.nmol)):
            nsites_total += self.dice.nmol[i] * len(self.system.molecule[i].atom)

        thickness = []
        picked_mols = []

        for proc in range(1, self.nprocs + 1):  # Run over folders

            simdir = "simfiles"
            path = (
                simdir
                + os.sep
                + "step{:02d}".format(cycle)
                + os.sep
                + "p{:02d}".format(proc)
            )
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
                    sys.exit("Error: wrong number of sites in file {}".format(file))

                box = xyzfile.pop(0).split()[-3:]
                box = [float(box[0]), float(box[1]), float(box[2])]
                sizes = self.system.molecule[0].sizes_of_molecule()
                thickness.append(
                    min(
                        [
                            (box[0] - sizes[0]) / 2,
                            (box[1] - sizes[1]) / 2,
                            (box[2] - sizes[2]) / 2,
                        ]
                    )
                )

                # Skip the first (reference) molecule
                xyzfile = xyzfile[nsitesref:]
                mol_count = 0
                for type in range(len(self.dice.nmol)):  # Run over types of molecules

                    if type == 0:
                        nmols = self.dice.nmol[0] - 1
                    else:
                        nmols = self.dice.nmol[type]

                    for mol in range(nmols):  # Run over molecules of each type

                        new_molecule = Molecule(self.system.molecule[type].molname)
                        # Run over sites of each molecule
                        for site in range(len(self.system.molecule[types].atom)):

                            # new_molecule.append({})
                            line = xyzfile.pop(0).split()

                            if (
                                line[0].title()
                                != atomsymb[
                                    self.system.molecule[type].atom[site].na.strip()
                                ]
                            ):
                                sys.exit("Error reading file {}".format(file))

                            new_molecule.add_atom(
                                Atom(
                                    self.system.molecule[type].atom[site].lbl,
                                    self.system.molecule[type].atom[site].na,
                                    self.system.molecule[type]
                                    .atom[site]
                                    .float(line[1]),
                                    self.system.molecule[type]
                                    .atom[site]
                                    .float(line[2]),
                                    self.system.molecule[type]
                                    .atom[site]
                                    .float(line[3]),
                                    self.system.molecule[type].atom[site].chg,
                                    self.system.molecule[type].atom[site].eps,
                                    self.system.molecule[type].atom[site].sig,
                                )
                            )

                        dist = self.system.molecule[0].minimum_distance(new_molecule)
                        if dist < thickness[-1]:
                            mol_count += 1
                            for atom in new_molecule.atom:
                                asec_charges.append({})
                                # vdw_meanfield.append({})

                                asec_charges[-1]["rx"] = atom.rx
                                asec_charges[-1]["ry"] = atom.ry
                                asec_charges[-1]["rz"] = atom.rz
                                asec_charges[-1]["chg"] = atom.chg / norm_factor

                                # if self.vdwforces == "yes":
                                #     vdw_meanfield[-1]["rx"] = atom["rx"]
                                #     vdw_meanfield[-1]["ry"] = atom["ry"]
                                #     vdw_meanfield[-1]["rz"] = atom["rz"]
                                #     vdw_meanfield[-1]["eps"] = atom["eps"]
                                #     vdw_meanfield[-1]["sig"] = atom["sig"]

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
            sum(picked_mols) / norm_factor
        )
        string += "were selected from each of the {} configurations ".format(
            len(picked_mols)
        )
        string += (
            "of the production simulations to form the ASEC, comprising a shell with "
        )
        string += "minimum thickness of {:>6.2f} Angstrom\n".format(
            sum(thickness) / norm_factor
        )

        self.outfile.write(textwrap.fill(string, 86))
        self.outfile.write("\n")

        otherfh = open("ASEC.dat", "w")
        for charge in asec_charges:
            otherfh.write(
                "{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
                    charge["rx"], charge["ry"], charge["rz"], charge["chg"]
                )
            )
        otherfh.close()

        return asec_charges