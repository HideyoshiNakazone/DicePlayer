from diceplayer.config.dice_config import DiceConfig
from diceplayer.config.gaussian_config import GaussianConfig
from diceplayer.config.player_config import PlayerConfig

import unittest


def get_config_dict():
    return {
        "opt": True,
        "mem": 12,
        "maxcyc": 100,
        "nprocs": 4,
        "ncores": 4,
        "dice": {
            "ljname": "test",
            "outname": "test",
            "dens": 1.0,
            "nmol": [1],
            "nstep": [1, 1],
        },
        "gaussian": {
            "level": "test",
            "qmprog": "g16",
            "keywords": "test",
        },
    }


class TestPlayerDTO(unittest.TestCase):
    def setUp(self) -> None:
        self.dice_dto = DiceConfig(
            ljname="test",
            outname="test",
            dens=1.0,
            nmol=[1],
            nstep=[1, 1],
        )
        self.gaussian_dto = GaussianConfig(
            level="test",
            qmprog="g16",
            keywords="test",
        )

    def test_class_instantiation(self):
        player_dto = PlayerConfig(
            opt=True,
            mem=12,
            maxcyc=100,
            nprocs=4,
            ncores=4,
            dice=self.dice_dto,
            gaussian=self.gaussian_dto,
        )

        self.assertIsInstance(player_dto, PlayerConfig)
        self.assertIsInstance(player_dto.dice, DiceConfig)
        self.assertIsInstance(player_dto.gaussian, GaussianConfig)

    def test_min_altsteps(self):
        player_dto = PlayerConfig(
            opt=True,
            mem=12,
            maxcyc=100,
            nprocs=4,
            ncores=4,
            altsteps=100,
            dice=self.dice_dto,
            gaussian=self.gaussian_dto,
        )

        self.assertEqual(player_dto.altsteps, 20000)

    def test_from_dict(self):
        player_dto = PlayerConfig.from_dict(get_config_dict())

        self.assertIsInstance(player_dto, PlayerConfig)
        self.assertIsInstance(player_dto.dice, DiceConfig)
        self.assertIsInstance(player_dto.gaussian, GaussianConfig)


if __name__ == "__main__":
    unittest.main()
