from diceplayer.shared.utils.dataclass_protocol import Dataclass

from abc import ABC, abstractmethod


class External(ABC):
    __slots__ = [
        'config'
    ]

    @abstractmethod
    def __init__(self, data: dict):
        pass

    @staticmethod
    @abstractmethod
    def set_config(data: dict) -> Dataclass:
        pass

    @abstractmethod
    def start(self, cycle: int):
        pass

    @abstractmethod
    def reset(self):
        pass
