from diceplayer.shared.config.player_dto import PlayerDTO

import unittest


class TestPlayerDTO(unittest.TestCase):
    def test_class_instantiation(self):
        player_dto = PlayerDTO(opt=True, maxcyc=100, nprocs=4)

        self.assertIsInstance(player_dto, PlayerDTO)

    def test_min_altsteps(self):
        player_dto = PlayerDTO(opt=True, maxcyc=100, nprocs=4, altsteps=100)

        self.assertEqual(player_dto.altsteps, 20000)

    def test_from_dict(self):
        player_dto = PlayerDTO.from_dict({
            'opt': True,
            'maxcyc': 100,
            'nprocs': 4
        })

        self.assertIsInstance(player_dto, PlayerDTO)


if __name__ == '__main__':
    unittest.main()
