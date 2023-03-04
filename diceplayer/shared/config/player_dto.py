from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dataclasses import dataclass
from dacite import from_dict


@dataclass
class PlayerDTO(Dataclass):
    opt: bool
    maxcyc: int
    nprocs: int

    qmprog: str = 'g16'
    altsteps: int = 20000
    simulation_dir = 'simfiles'

    def __post_init__(self):
        MIN_STEP = 20000
        # altsteps value is always the nearest multiple of 1000
        self.altsteps = round(max(MIN_STEP, self.altsteps) / 1000) * 1000

    @classmethod
    def from_dict(cls, param: dict):
        return from_dict(PlayerDTO, param)
