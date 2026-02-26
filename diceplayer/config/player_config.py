from diceplayer.config.dice_config import DiceConfig
from diceplayer.config.gaussian_config import GaussianConfig

from pydantic import BaseModel, Field, model_validator
from typing_extensions import Self

from pathlib import Path


MIN_STEP = 20000


class PlayerConfig(BaseModel):
    """
    Data Transfer Object for the player configuration.
    """

    opt: bool = Field(..., description="Whether to perform geometry optimization")
    maxcyc: int = Field(
        ..., description="Maximum number of cycles for the geometry optimization"
    )
    nprocs: int = Field(
        ..., description="Number of processors to use for the QM calculations"
    )
    ncores: int = Field(
        ..., description="Number of cores to use for the QM calculations"
    )

    dice: DiceConfig = Field(..., description="Dice configuration")
    gaussian: GaussianConfig = Field(..., description="Gaussian configuration")

    mem: int = Field(None, description="Memory configuration")
    switchcyc: int = Field(3, description="Switch cycle configuration")
    qmprog: str = Field("g16", description="QM program to use for the calculations")
    altsteps: int = Field(
        20000, description="Number of steps for the alternate simulation"
    )
    geoms_file: Path = Field(
        "geoms.xyz", description="File name for the geometries output"
    )
    simulation_dir: Path = Field(
        "simfiles", description="Directory name for the simulation files"
    )

    @model_validator(mode="after")
    def validate_altsteps(self) -> Self:
        self.altsteps = round(max(MIN_STEP, self.altsteps) / 1000) * 1000
        return self
