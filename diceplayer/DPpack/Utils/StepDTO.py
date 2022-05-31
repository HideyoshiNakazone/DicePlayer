from dataclasses import dataclass
from typing import List

from diceplayer.DPpack.Environment.Molecule import Molecule


@dataclass
class StepDTO:

    initcyc: int
    nprocs: int
    ncores: int
    altsteps: int
    switchcyc: int
    opt: str
    nmol: List[int]
    molecule: List[Molecule]