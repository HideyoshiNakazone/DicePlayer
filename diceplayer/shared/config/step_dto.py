from diceplayer.shared.environment.molecule import Molecule

from dataclasses import dataclass
from typing import List


@dataclass
class StepDTO:
    nprocs: int = None
    ncores: int = None
    altsteps: int = None
    switchcyc: int = None
    opt: str = None
    nmol: List[int] = None
    molecule: List[Molecule] = None

    charges: List[float] = None
    position: List[float] = None