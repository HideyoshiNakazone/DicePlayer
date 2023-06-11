from diceplayer.shared.config.dice_config import DiceConfig
from diceplayer.shared.config.gaussian_config import GaussianDTO
from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dacite import from_dict

from dataclasses import dataclass


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
    gaussian: GaussianDTO

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
    def from_dict(cls, param: dict):
        if param["dice"] is None:
            raise ValueError("Error: 'dice' keyword not specified in config file.")
        param["dice"] = DiceConfig.from_dict(param["dice"])

        if param["gaussian"] is None:
            raise ValueError("Error: 'gaussian' keyword not specified in config file.")
        param["gaussian"] = GaussianDTO.from_dict(param["gaussian"])

        return from_dict(PlayerConfig, param)
