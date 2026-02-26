from typing import Literal

from pydantic import BaseModel, Field

from diceplayer.shared.utils.dataclass_protocol import Dataclass

from dataclasses import dataclass, fields


class GaussianConfig(BaseModel):
    """
    Data Transfer Object for the Gaussian configuration.
    """

    level: str = Field(..., description="Level of theory for the QM calculations")
    qmprog: Literal["g03", "g09", "g16"] = Field("g16", description="QM program to use for the calculations")

    chgmult: list[int] = Field(default_factory=lambda: [0, 1], description="List of charge and multiplicity for the QM calculations")
    pop: str = Field("chelpg", description="Population analysis method for the QM calculations")
    chg_tol: float = Field(0.01, description="Charge tolerance for the QM calculations")
    keywords: str = Field(None, description="Additional keywords for the QM calculations")
