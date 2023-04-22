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

        if not isinstance(self.ljname, str):
            raise ValueError(
                "Error: 'ljname' keyword not specified in config file"
            )

        if not isinstance(self.outname, str):
            raise ValueError(
                "Error: 'outname' keyword not specified in config file"
            )

        if not isinstance(self.dens, float):
            raise ValueError(
                "Error: 'dens' keyword not specified in config file"
            )

        if not isinstance(self.nmol, list):
            raise ValueError(
                "Error: 'nmol' keyword not defined appropriately in config file"
            )

        if not isinstance(self.nstep, list):
            raise ValueError(
                "Error: 'nstep' keyword not defined appropriately in config file"
            )

    @classmethod
    def from_dict(cls, param: dict):
        return from_dict(DiceDTO, param)
