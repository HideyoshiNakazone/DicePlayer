from diceplayer import logger
from diceplayer.shared.config.dice_config import DiceConfig
from diceplayer.shared.config.gaussian_config import GaussianDTO
from diceplayer.shared.config.player_config import PlayerConfig
from diceplayer.shared.environment.atom import Atom
from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.environment.system import System
from diceplayer.shared.interface.dice_interface import DiceInterface
from diceplayer.shared.interface.gaussian_interface import GaussianInterface
from diceplayer.shared.utils.dataclass_protocol import Dataclass
from diceplayer.shared.utils.misc import weekday_date_time
from diceplayer.shared.utils.ptable import atomsymb

import yaml

import os
import pickle
import sys
from dataclasses import fields
from pathlib import Path
from typing import Tuple, Type

ENV = ["OMP_STACKSIZE"]


class Player:
    def __init__(self, infile: str = None, optimization: bool = False):
        if infile is None and optimization is False:
            raise ValueError("Must specify either infile or optimization")

        elif infile is not None:
            self.config = self.set_config(self.read_keywords(infile))

            self.system = System()

            self.initial_cycle = 1

        elif optimization is True:
            save = self.load_run_from_pickle()

            self.config = save[0]

            self.system = save[1]

            self.initial_cycle = save[2] + 1

        else:
            raise ValueError("Must specify either infile or config")

        self.dice_interface = DiceInterface()
        self.gaussian_interface = GaussianInterface()

    def start(self):
        logger.info(
            "==========================================================================================\n"
            "Starting the iterative process.\n"
            "==========================================================================================\n"
        )

        for cycle in range(self.initial_cycle, self.initial_cycle + self.config.maxcyc):
            logger.info(
                f"------------------------------------------------------------------------------------------\n"
                f"                                         Step # {cycle}\n"
                f"------------------------------------------------------------------------------------------\n"
            )

            self.dice_start(cycle)

            try:
                self.gaussian_start(cycle)
            except StopIteration:
                break

            self.save_run_in_pickle(cycle)

    def prepare_system(self):
        for i, mol in enumerate(self.system.molecule):
            logger.info(f"Molecule {i + 1} - {mol.molname}")

            mol.print_mol_info()
            logger.info(
                "\n    Translating and rotating molecule to standard orientation..."
            )

            mol.standard_orientation()
            logger.info("\n Done")
            logger.info("\nNew values:\n")
            mol.print_mol_info()

            logger.info("\n")

    def create_simulation_dir(self):
        simulation_dir_path = Path(self.config.simulation_dir)
        if simulation_dir_path.exists():
            raise FileExistsError(
                f"Error: a file or a directory {self.config.simulation_dir} already exists,"
                f" move or delete the simfiles directory to continue."
            )
        simulation_dir_path.mkdir()

    def create_geoms_file(self):
        geoms_file_path = Path(self.config.geoms_file)
        if geoms_file_path.exists():
            raise FileExistsError(
                f"Error: a file or a directory {self.config.geoms_file} already exists,"
                f" move or delete the simfiles directory to continue."
            )
        geoms_file_path.touch()

    def print_keywords(self) -> None:
        def log_keywords(config: Dataclass, dto: Type[Dataclass]):
            for key in sorted(list(map(lambda f: f.name, fields(dto)))):
                if getattr(config, key) is not None:
                    if isinstance(getattr(config, key), list):
                        string = " ".join(str(x) for x in getattr(config, key))
                        logger.info(f"{key} = [ {string} ]")
                    else:
                        logger.info(f"{key} = {getattr(config, key)}")

        logger.info(
            "##########################################################################################\n"
            "#############               Welcome to DICEPLAYER version 1.0                #############\n"
            "##########################################################################################\n"
        )
        logger.info("Your python version is {}\n".format(sys.version))
        logger.info("Program started on {}\n".format(weekday_date_time()))
        logger.info("Environment variables:")
        for var in ENV:
            logger.info(
                "{} = {}\n".format(
                    var, (os.environ[var] if var in os.environ else "Not set")
                )
            )

        logger.info(
            "------------------------------------------------------------------------------------------\n"
            "                         DICE variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
        )

        log_keywords(self.config.dice, DiceConfig)

        logger.info(
            "------------------------------------------------------------------------------------------\n"
            "                         GAUSSIAN variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
        )

        log_keywords(self.config.gaussian, GaussianDTO)

        logger.info("\n")

    def read_potentials(self):
        ljname_path = Path(self.config.dice.ljname)
        if ljname_path.exists():
            with open(self.config.dice.ljname) as file:
                ljc_data = file.readlines()
        else:
            raise RuntimeError(f"Potential file {self.config.dice.ljname} not found.")

        combrule = ljc_data.pop(0).split()[0]
        if combrule not in ("*", "+"):
            sys.exit(
                "Error: expected a '*' or a '+' sign in 1st line of file {}".format(
                    self.config.dice.ljname
                )
            )
        self.config.dice.combrule = combrule

        ntypes = ljc_data.pop(0).split()[0]
        if not ntypes.isdigit():
            sys.exit(
                "Error: expected an integer in the 2nd line of file {}".format(
                    self.config.dice.ljname
                )
            )
        ntypes = int(ntypes)

        if ntypes != len(self.config.dice.nmol):
            sys.exit(
                f"Error: number of molecule types in file {self.config.dice.ljname} "
                f"must match that of 'nmol' keyword in config file"
            )

        for i in range(ntypes):
            try:
                nsites, molname = ljc_data.pop(0).split()[:2]
            except ValueError:
                raise ValueError(
                    f"Error: expected nsites and molname for the molecule type {i + 1}"
                )

            if not nsites.isdigit():
                raise ValueError(
                    f"Error: expected nsites to be an integer for molecule type {i + 1}"
                )

            nsites = int(nsites)
            self.system.add_type(Molecule(molname))

            atom_fields = ["lbl", "na", "rx", "ry", "rz", "chg", "eps", "sig"]
            for j in range(nsites):
                new_atom = dict(zip(atom_fields, ljc_data.pop(0).split()))
                self.system.molecule[i].add_atom(
                    Atom(**self.validate_atom_dict(i, j, new_atom))
                )

    def print_potentials(self) -> None:
        formatstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f} {:>9.4f}"
        logger.info(
            "==========================================================================================\n"
            f"                    Potential parameters from file {self.config.dice.ljname}:\n"
            "------------------------------------------------------------------------------------------"
            "\n"
        )

        logger.info(f"Combination rule: {self.config.dice.combrule}")
        logger.info(f"Types of molecules: {len(self.system.molecule)}\n")

        i = 0
        for mol in self.system.molecule:
            i += 1
            logger.info("{} atoms in molecule type {}:".format(len(mol.atom), i))
            logger.info(
                "---------------------------------------------------------------------------------"
            )
            logger.info(
                "Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass"
            )
            logger.info(
                "---------------------------------------------------------------------------------"
            )

            for atom in mol.atom:
                logger.info(
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

            logger.info("\n")

    def dice_start(self, cycle: int):
        self.dice_interface.configure(
            self.config,
            self.system,
        )

        self.dice_interface.start(cycle)

        self.dice_interface.reset()

    def gaussian_start(self, cycle: int):
        self.gaussian_interface.configure(
            self.config,
            self.system,
        )

        result = self.gaussian_interface.start(cycle)

        self.gaussian_interface.reset()

        if self.config.opt:
            if "position" not in result:
                raise RuntimeError("Optimization failed. No position found in result.")

            self.system.update_molecule(result["position"])

        else:
            if "charges" not in result:
                raise RuntimeError(
                    "Charges optimization failed. No charges found in result."
                )

            diff = self.system.molecule[0].update_charges(result["charges"])

            self.system.print_charges_and_dipole(cycle)
            self.print_geoms(cycle)

            if diff < self.config.gaussian.chg_tol:
                logger.info(f"Charges converged after {cycle} cycles.")
                raise StopIteration()

    def print_geoms(self, cycle: int):
        with open(self.config.geoms_file, "a") as file:
            file.write(f"Cycle # {cycle}\n")

            for atom in self.system.molecule[0].atom:
                symbol = atomsymb[atom.na]
                file.write(
                    f"{symbol:<2s}    {atom.rx:>10.6f}  {atom.ry:>10.6f}  {atom.rz:>10.6f}\n"
                )

            file.write("\n")

    @staticmethod
    def validate_atom_dict(molecule_type, molecule_site, atom_dict: dict) -> dict:
        molecule_type += 1
        molecule_site += 1

        if len(atom_dict) < 8:
            raise ValueError(
                f"Invalid number of fields for site {molecule_site} for molecule type {molecule_type}."
            )

        try:
            atom_dict["lbl"] = int(atom_dict["lbl"])
        except Exception:
            raise ValueError(
                f"Invalid lbl fields for site {molecule_site} for molecule type {molecule_type}."
            )

        try:
            atom_dict["na"] = int(atom_dict["na"])
        except Exception:
            raise ValueError(
                f"Invalid na fields for site {molecule_site} for molecule type {molecule_type}."
            )

        try:
            atom_dict["rx"] = float(atom_dict["rx"])
        except Exception:
            raise ValueError(
                f"Invalid rx fields for site {molecule_site} for molecule type {molecule_type}. "
                f"Value must be a float."
            )

        try:
            atom_dict["ry"] = float(atom_dict["ry"])
        except Exception:
            raise ValueError(
                f"Invalid ry fields for site {molecule_site} for molecule type {molecule_type}. "
                f"Value must be a float."
            )

        try:
            atom_dict["rz"] = float(atom_dict["rz"])
        except Exception:
            raise ValueError(
                f"Invalid rz fields for site {molecule_site} for molecule type {molecule_type}. "
                f"Value must be a float."
            )

        try:
            atom_dict["chg"] = float(atom_dict["chg"])
        except Exception:
            raise ValueError(
                f"Invalid chg fields for site {molecule_site} for molecule type {molecule_type}. "
                f"Value must be a float."
            )

        try:
            atom_dict["eps"] = float(atom_dict["eps"])
        except Exception:
            raise ValueError(
                f"Invalid eps fields for site {molecule_site} for molecule type {molecule_type}. "
                f"Value must be a float."
            )

        try:
            atom_dict["sig"] = float(atom_dict["sig"])
        except Exception:
            raise ValueError(
                f"Invalid sig fields for site {molecule_site} for molecule type {molecule_type}. "
                f"Value must be a float."
            )

        return atom_dict

    def print_results(self):
        formatstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f} {:>9.4f}"

        mol = self.system.molecule[0]
        logger.info("{} atoms in molecule type {}:".format(len(mol.atom), 1))
        logger.info(
            "---------------------------------------------------------------------------------"
        )
        logger.info(
            "Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass"
        )
        logger.info(
            "---------------------------------------------------------------------------------"
        )

        for atom in mol.atom:
            logger.info(
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

        logger.info("\n")

    def save_run_in_pickle(self, cycle):
        try:
            with open("latest-step.pkl", "wb") as pickle_file:
                pickle.dump((self.config, self.system, cycle), pickle_file)
        except Exception:
            raise RuntimeError(f"Could not save pickle file latest-step.pkl.")

    @staticmethod
    def load_run_from_pickle() -> Tuple[PlayerConfig, System, int]:
        pickle_path = Path("latest-step.pkl")
        try:
            with open(pickle_path, "rb") as pickle_file:
                save = pickle.load(pickle_file)
            return save[0], save[1], save[2] + 1

        except Exception:
            raise RuntimeError(f"Could not load pickle file {pickle_path}.")

    @staticmethod
    def set_config(data: dict) -> PlayerConfig:
        return PlayerConfig.from_dict(data)

    @staticmethod
    def read_keywords(infile) -> dict:
        with open(infile, "r") as yml_file:
            config = yaml.load(yml_file, Loader=yaml.SafeLoader)

        if "diceplayer" in config:
            return config.get("diceplayer")

        raise RuntimeError(f"Could not find diceplayer section in {infile}.")

    @classmethod
    def from_file(cls, infile: str) -> "Player":
        return cls(infile=infile)

    @classmethod
    def from_save(cls):
        return cls(optimization=True)
