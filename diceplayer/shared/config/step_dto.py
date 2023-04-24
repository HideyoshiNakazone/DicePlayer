from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.config.player_dto import PlayerDTO

from dataclasses import dataclass
from typing import List


@dataclass
class StepDTO:
    """
    Data Transfer Object for the step configuration.
    """
    ncores: int
    nprocs: int
    simulation_dir: str

    altsteps: int

    nmol: List[int] = None
    molecule: List[Molecule] = None
    charges: List[float] = None
    position: List[float] = None
