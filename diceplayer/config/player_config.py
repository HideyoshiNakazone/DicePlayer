from diceplayer.config.dice_config import DiceConfig
from diceplayer.config.gaussian_config import GaussianConfig
from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dataclasses import dataclass, fields


@dataclass
class PlayerConfig(Dataclass):
    """
    Data Transfer Object for the player configuration.
    """

    opt: bool
    maxcyc: int
    nprocs: int
    ncores: int

    dice: DiceConfig
    gaussian: GaussianConfig

    mem: int = None
    switchcyc: int = 3
    qmprog: str = "g16"
    altsteps: int = 20000
    geoms_file = "geoms.xyz"
    simulation_dir = "simfiles"

    def __post_init__(self):
        MIN_STEP = 20000
        # altsteps value is always the nearest multiple of 1000
        self.altsteps = round(max(MIN_STEP, self.altsteps) / 1000) * 1000

    @classmethod
    def from_dict(cls, params: dict):
        if params["dice"] is None:
            raise ValueError("Error: 'dice' keyword not specified in config file.")
        params["dice"] = DiceConfig.model_validate(params["dice"])

        if params["gaussian"] is None:
            raise ValueError("Error: 'gaussian' keyword not specified in config file.")
        params["gaussian"] = GaussianConfig.model_validate(params["gaussian"])

        params = {f.name: params[f.name] for f in fields(cls) if f.name in params}

        return cls(**params)
