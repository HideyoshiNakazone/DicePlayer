from diceplayer.shared.config.dice_dto import DiceDTO

import unittest


class TestDiceDto(unittest.TestCase):
    def test_class_instantiation(self):
        dice_dto = DiceDTO(
            ljname='test',
            outname='test',
            ncores=1,
            dens=1.0,
            nmol=[1],
            nstep=[1],
        )

        self.assertIsInstance(dice_dto, DiceDTO)

    def test_validate_jname(self):
        with self.assertRaises(ValueError) as ex:
            DiceDTO(
                ljname=None,
                outname='test',
                ncores=1,
                dens=1.0,
                nmol=[1],
                nstep=[1],
            )
            self.assertEqual(ex.exception, "Error: 'ljname' keyword not specified in config file")

    def test_validate_outname(self):
        with self.assertRaises(ValueError) as ex:
            DiceDTO(
                ljname='test',
                outname=None,
                ncores=1,
                dens=1.0,
                nmol=[1],
                nstep=[1],
            )
            self.assertEqual(ex.exception, "Error: 'outname' keyword not specified in config file")

    def test_validate_dens(self):
        with self.assertRaises(ValueError) as ex:
            DiceDTO(
                ljname='test',
                outname='test',
                ncores=1,
                dens=None,
                nmol=[1],
                nstep=[1],
            )
            self.assertEqual(ex.exception, "Error: 'dens' keyword not specified in config file")

    def test_validate_nmol(self):
        with self.assertRaises(ValueError) as ex:
            DiceDTO(
                ljname='test',
                outname='test',
                ncores=1,
                dens=1.0,
                nmol=0,
                nstep=[1],
            )
            self.assertEqual(ex.exception, "Error: 'nmol' keyword not defined appropriately in config file")

    def test_validate_nstep(self):
        with self.assertRaises(ValueError) as ex:
            DiceDTO(
                ljname='test',
                outname='test',
                ncores=1,
                dens=1.0,
                nmol=[1],
                nstep=0,
            )
            self.assertEqual(ex.exception, "Error: 'nstep' keyword not defined appropriately in config file")

    def test_from_dict(self):
        dice_dto = DiceDTO.from_dict({
            'ljname': 'test',
            'outname': 'test',
            'ncores': 1,
            'dens': 1.0,
            'nmol': [1],
            'nstep': [1],
        })

        self.assertIsInstance(dice_dto, DiceDTO)