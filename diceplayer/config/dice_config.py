from diceplayer.shared.utils.dataclass_protocol import Dataclass

from pydantic import BaseModel, Field
from typing_extensions import List, Literal

from dataclasses import dataclass, fields


class DiceConfig(BaseModel):
    """
    Data Transfer Object for the Dice configuration.
    """

    ljname: str = Field(..., description="Name of the Lennard-Jones potential file")
    outname: str = Field(
        ..., description="Name of the output file for the simulation results"
    )
    dens: float = Field(..., description="Density of the system")
    nmol: List[int] = Field(
        ..., description="List of the number of molecules for each component"
    )
    nstep: List[int] = Field(
        ...,
        description="List of the number of steps for each component",
        min_length=2,
        max_length=3,
    )

    upbuf: int = Field(
        360, description="Buffer size for the potential energy calculations"
    )
    combrule: Literal["+", "*"] = Field(
        "*", description="Combination rule for the Lennard-Jones potential"
    )
    isave: int = Field(1000, description="Frequency of saving the simulation results")
    press: float = Field(1.0, description="Pressure of the system")
    temp: float = Field(300.0, description="Temperature of the system")
    progname: str = Field(
        "dice", description="Name of the program to run the simulation"
    )
    randominit: str = Field(
        "first", description="Method for initializing the random number generator"
    )
    seed: int | None = Field(None, description="Seed for the random number generator")
