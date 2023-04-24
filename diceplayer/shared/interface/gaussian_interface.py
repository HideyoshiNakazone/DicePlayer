from diceplayer.shared.config.gaussian_dto import GaussianDTO
from diceplayer.shared.interface import Interface


class GaussianInterface(Interface):

    def __init__(self, data: dict):
        self.config: GaussianDTO = self.set_config(data)

    @staticmethod
    def set_config(data: dict) -> GaussianDTO:
        return GaussianDTO.from_dict(data)

    def configure(self):
        pass

    def start(self):
        pass

    def reset(self):
        pass

