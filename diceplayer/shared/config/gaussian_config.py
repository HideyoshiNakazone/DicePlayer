from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dataclasses import dataclass, fields


@dataclass
class GaussianDTO(Dataclass):
    """
    Data Transfer Object for the Gaussian configuration.
    """

    level: str
    qmprog: str

    chgmult = [0, 1]
    pop: str = "chelpg"
    chg_tol: float = 0.01
    keywords: str = None

    def __post_init__(self):
        if self.qmprog not in ("g03", "g09", "g16"):
            raise ValueError("Error: invalid qmprog value.")
        if self.level is None:
            raise ValueError("Error: 'level' keyword not specified in config file.")

    @classmethod
    def from_dict(cls, params: dict):
        params = {f.name: params[f.name] for f in fields(cls) if f.name in params}

        return cls(**params)
