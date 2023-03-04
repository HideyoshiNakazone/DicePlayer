from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dataclasses import dataclass
from dacite import from_dict


@dataclass
class GaussianDTO(Dataclass):
    level: str
    qmprog: str
    keywords: str

    chgmult = [0, 1]
    pop: str = 'chelpg'

    def __post_init__(self):
        if self.qmprog not in ("g03", "g09", "g16"):
            raise ValueError(
                "Error: invalid qmprog value."
            )
        if self.level is None:
            raise ValueError(
                "Error: 'level' keyword not specified in config file."
            )

    @classmethod
    def from_dict(cls, param: dict):
        return from_dict(GaussianDTO, param)
