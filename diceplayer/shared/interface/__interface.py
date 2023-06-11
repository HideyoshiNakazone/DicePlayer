from __future__ import annotations

from diceplayer.shared.config.player_config import PlayerConfig
from diceplayer.shared.environment.system import System

from abc import ABC, abstractmethod


class Interface(ABC):
    __slots__ = ["step", "system"]

    def __init__(self):
        self.system: System | None = None
        self.step: PlayerConfig | None = None

    @abstractmethod
    def configure(self, step: PlayerConfig, system: System):
        pass

    @abstractmethod
    def start(self, cycle: int):
        pass

    @abstractmethod
    def reset(self):
        pass
