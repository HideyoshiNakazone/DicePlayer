from diceplayer.shared.utils.dataclass_protocol import Dataclass
from diceplayer.shared.external.__external import External
from diceplayer.shared.config.dice_dto import DiceDTO


class Dice(External):

    def __init__(self, data: dict):
        self.config: DiceDTO = self.set_config(data)

    @staticmethod
    def set_config(data: dict) -> DiceDTO:
        return DiceDTO.from_dict(data)

    def configure(self):
        pass

    def start(self):
        pass

    def reset(self):
        pass
