from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dacite import from_dict

from dataclasses import dataclass
from typing import List


@dataclass
class DiceConfig(Dataclass):
    """
    Data Transfer Object for the Dice configuration.
    """

    ljname: str
    outname: str
    dens: float
    nmol: List[int]
    nstep: List[int]

    upbuf = 360
    combrule = "*"
    isave: int = 1000
    press: float = 1.0
    temp: float = 300.0
    progname: str = "dice"
    randominit: str = "first"

    def __post_init__(self):
        if not isinstance(self.ljname, str):
            raise ValueError("Error: 'ljname' keyword not specified in config file")

        if not isinstance(self.outname, str):
            raise ValueError("Error: 'outname' keyword not specified in config file")

        if not isinstance(self.dens, float):
            raise ValueError("Error: 'dens' keyword not specified in config file")

        if not isinstance(self.nmol, list):
            raise ValueError(
                "Error: 'nmol' keyword not defined appropriately in config file"
            )

        if not isinstance(self.nstep, list) or len(self.nstep) not in (2, 3):
            raise ValueError(
                "Error: 'nstep' keyword not defined appropriately in config file"
            )

    @classmethod
    def from_dict(cls, param: dict):
        return from_dict(DiceConfig, param)
