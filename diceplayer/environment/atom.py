from diceplayer.utils.ptable import PTable, AtomInfo

from dataclasses import dataclass


@dataclass
class Atom:
    """
    Atom class declaration. This class is used throughout the DicePlayer program to represent atoms.
    """

    lbl: int
    na: int
    rx: float
    ry: float
    rz: float
    chg: float
    eps: float
    sig: float

    @property
    def mass(self) -> float:
        return PTable.get_atomic_mass(self.na)

    @property
    def atom_info(self) -> AtomInfo:
        return PTable.get_from_atomic_number(self.na)
