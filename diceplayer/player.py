from diceplayer.shared.interface.gaussian_interface import GaussianInterface
from diceplayer.shared.interface.dice_interface import DiceInterface
from diceplayer.shared.utils.dataclass_protocol import Dataclass
from diceplayer.shared.config.gaussian_dto import GaussianDTO
from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.utils.misc import weekday_date_time
from diceplayer.shared.config.player_dto import PlayerDTO
from diceplayer.shared.environment.system import System
from diceplayer.shared.config.step_dto import StepDTO
from diceplayer.shared.config.dice_dto import DiceDTO
from diceplayer.shared.environment.atom import Atom
from diceplayer.shared.utils.ptable import atommass

from dataclasses import fields
from pathlib import Path
from typing import Type
import logging
import yaml
import sys
import os


ENV = ["OMP_STACKSIZE"]


class Player:
    __slots__ = [
        'config',
        'system',
        'dice',
        'gaussian',
    ]

    def __init__(self, infile: str):
        config_data = self.read_keywords(infile)

        self.system = System()

        self.config = self.set_config(
            config_data.get("diceplayer")
        )

        self.gaussian = GaussianInterface(config_data.get("gaussian"))
        self.dice = DiceInterface(config_data.get("dice"))

    def start(self):
        self.print_keywords()

        self.create_simulation_dir()

        self.read_potentials()
        # self.print_potentials()

        self.dice_start(1)
        self.dice_start(2)

    def create_simulation_dir(self):
        simulation_dir_path = Path(self.config.simulation_dir)
        if simulation_dir_path.exists():
            raise FileExistsError(
                f"Error: a file or a directory {self.config.simulation_dir} already exists,"
                f" move or delete the simfiles directory to continue."
            )
        try:
            simulation_dir_path.mkdir()
        except FileExistsError:
            OSError(
                f"Error: cannot make directory {self.config.simulation_dir}"
            )

    def print_keywords(self) -> None:

        def log_keywords(config: Dataclass, dto: Type[Dataclass]):
            for key in sorted(list(map(lambda f: f.name, fields(dto)))):
                if getattr(config, key) is not None:
                    if isinstance(getattr(config, key), list):
                        string = " ".join(str(x) for x in getattr(config, key))
                        logging.info(f"{key} = [ {string} ]")
                    else:
                        logging.info(f"{key} = {getattr(config, key)}")

        logging.info(
            "##########################################################################################\n"
            "#############               Welcome to DICEPLAYER version 1.0                #############\n"
            "##########################################################################################\n"
            "\n"
        )
        logging.info("Your python version is {}\n".format(sys.version))
        logging.info("\n")
        logging.info("Program started on {}\n".format(weekday_date_time()))
        logging.info("\n")
        logging.info("Environment variables:\n")
        for var in ENV:
            logging.info(
                "{} = {}\n".format(
                    var, (os.environ[var] if var in os.environ else "Not set")
                )
            )

        logging.info(
            "\n==========================================================================================\n"
            "                         CONTROL variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
            "\n"
        )

        logging.info("\n")

        logging.info(
            "------------------------------------------------------------------------------------------\n"
            "                         DICE variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
            "\n"
        )

        log_keywords(self.dice.config, DiceDTO)

        logging.info("\n")

        logging.info(
            "------------------------------------------------------------------------------------------\n"
            "                         GAUSSIAN variables being used in this run:\n"
            "------------------------------------------------------------------------------------------\n"
            "\n"
        )

        log_keywords(self.gaussian.config, GaussianDTO)

        logging.info("\n")

    def read_potentials(self):
        try:
            with open(self.dice.config.ljname) as file:
                ljdata = file.readlines()
        except FileNotFoundError:
            raise RuntimeError(
                f"Potential file {self.dice.config.ljname} not found."
            )

        combrule = ljdata.pop(0).split()[0]
        if combrule not in ("*", "+"):
            sys.exit(
                "Error: expected a '*' or a '+' sign in 1st line of file {}".format(
                    self.dice.config.ljname
                )
            )
        self.dice.config.combrule = combrule

        ntypes = ljdata.pop(0).split()[0]
        if not ntypes.isdigit():
            sys.exit(
                "Error: expected an integer in the 2nd line of file {}".format(
                    self.dice.config.ljname
                )
            )
        ntypes = int(ntypes)

        if ntypes != len(self.dice.config.nmol):
            sys.exit(
                f"Error: number of molecule types in file {self.dice.config.ljname}"
                f"must match that of 'nmol' keyword in config file"
            )

        for i in range(ntypes):

            nsites, molname = ljdata.pop(0).split()[:2]

            if not nsites.isdigit():
                raise ValueError(
                    f"Error: expected nsites to be an integer for molecule type {i}"
                )

            if molname is None:
                raise ValueError(
                    f"Error: expected molecule name for molecule type {i}"
                )

            nsites = int(nsites)
            self.system.add_type(nsites, Molecule(molname))

            atom_fields = ["lbl", "na", "rx", "ry", "rz", "chg", "eps", "sig"]
            for j in range(nsites):
                new_atom = dict(zip(
                    atom_fields,
                    ljdata.pop(0).split()
                ))
                self.system.molecule[i].add_atom(
                    Atom(**self.validate_atom_dict(i, j, new_atom))
                )

    def dice_start(self, cycle: int):
        self.dice.configure(
            StepDTO(
                ncores=self.config.ncores,
                nprocs=self.config.nprocs,
                simulation_dir=self.config.simulation_dir,
                altsteps=self.config.altsteps,
                molecule=self.system.molecule,
                nmol=self.system.nmols,
            )
        )

        self.dice.start(cycle)

        self.dice.reset()

    def gaussian_start(self):
        self.gaussian.start()

    @staticmethod
    def validate_atom_dict(molecule_type, molecule_site, atom_dict: dict) -> dict:
        molecule_type += 1
        molecule_site += 1

        if len(atom_dict) < 8:
            raise ValueError(
                f'Invalid number of fields for site {molecule_site} for molecule type {molecule_type}.'
            )

        try:
            atom_dict['lbl'] = int(atom_dict['lbl'])
        except ValueError:
            raise ValueError(
                f'Invalid lbl fields for site {molecule_site} for molecule type {molecule_type}.'
            )

        try:
            atom_dict['na'] = int(atom_dict['na'])
        except ValueError:
            raise ValueError(
                f'Invalid na fields for site {molecule_site} for molecule type {molecule_type}.'
            )

        try:
            atom_dict['rx'] = float(atom_dict['rx'])
        except ValueError:
            raise ValueError(
                f'Invalid rx fields for site {molecule_site} for molecule type {molecule_type}.'
                f'Value must be a float.'
            )

        try:
            atom_dict['ry'] = float(atom_dict['ry'])
        except ValueError:
            raise ValueError(
                f'Invalid ry fields for site {molecule_site} for molecule type {molecule_type}.'
                f'Value must be a float.'
            )

        try:
            atom_dict['rz'] = float(atom_dict['rx'])
        except ValueError:
            raise ValueError(
                f'Invalid rz fields for site {molecule_site} for molecule type {molecule_type}.'
                f'Value must be a float.'
            )

        try:
            atom_dict['chg'] = float(atom_dict['chg'])
        except ValueError:
            raise ValueError(
                f'Invalid chg fields for site {molecule_site} for molecule type {molecule_type}.'
                f'Value must be a float.'
            )

        try:
            atom_dict['eps'] = float(atom_dict['eps'])
        except ValueError:
            raise ValueError(
                f'Invalid eps fields for site {molecule_site} for molecule type {molecule_type}.'
                f'Value must be a float.'
            )

        try:
            atom_dict['sig'] = float(atom_dict['sig'])
        except ValueError:
            raise ValueError(
                f'Invalid sig fields for site {molecule_site} for molecule type {molecule_type}.'
                f'Value must be a float.'
            )

        return atom_dict

    @staticmethod
    def set_config(data: dict) -> PlayerDTO:
        return PlayerDTO.from_dict(data)

    @staticmethod
    def read_keywords(infile) -> dict:
        with open(infile, 'r') as yml_file:
            return yaml.load(yml_file, Loader=yaml.SafeLoader)
