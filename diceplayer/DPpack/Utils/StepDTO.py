from dataclasses import dataclass
from typing import List

from diceplayer.DPpack.Environment.Molecule import Molecule


@dataclass
class StepDTO:

    cycle: int = None
    initcyc: int = None
    nprocs: int = None
    ncores: int = None
    altsteps: int = None
    switchcyc: int = None
    opt: str = None
    nmol: List[int] = None
    molecule: List[Molecule] = None

    charges: List[float] = None
    position: List[float] = None