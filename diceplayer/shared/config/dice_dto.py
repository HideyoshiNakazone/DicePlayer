from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dataclasses import dataclass
from dacite import from_dict
from typing import List


@dataclass
class DiceDTO(Dataclass):

    ljname: str
    outname: str
    ncores: int
    dens: float
    nmol: List[int]
    nstep: List[int]

    upbuf = 360
    combrule = "*"
    isave: int = 1000
    press: float = 1.0
    temp: float = 300.0
    randominit: str = 'first'

    def __post_init__(self):

        if self.ljname is None:
            raise ValueError(
                "Error: 'ljname' keyword not specified in config file"
            )

        if self.outname is None:
            raise ValueError(
                "Error: 'outname' keyword not specified in config file"
            )

        if self.dens is None:
            raise ValueError(
                "Error: 'dens' keyword not specified in config file"
            )

        if self.nmol == 0:
            raise ValueError(
                "Error: 'nmol' keyword not defined appropriately in config file"
            )

        if self.nstep == 0:
            raise ValueError(
                "Error: 'nstep' keyword not defined appropriately in config file"
            )

    @classmethod
    def from_dict(cls, param: dict):
        return from_dict(DiceDTO, param)
