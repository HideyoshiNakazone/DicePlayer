from diceplayer.shared.external.dice import Dice

from unittest import mock
import unittest


class TestDice(unittest.TestCase):
    def test_class_instantiation(self):
        dice = Dice(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1],
            }
        )

        self.assertIsInstance(dice, Dice)

    def test_configure(self):
        dice = Dice(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1],
            }
        )

        self.assertFalse(hasattr(dice, 'step'))

        dice.configure('test')

        self.assertTrue(hasattr(dice, 'step'))

    def test_reset(self):
        dice = Dice(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1],
            }
        )

        dice.configure('test')

        self.assertTrue(hasattr(dice, 'step'))

        dice.reset()

        self.assertFalse(hasattr(dice, 'step'))

    @mock.patch('diceplayer.shared.external.dice.connection')
    @mock.patch('diceplayer.shared.external.dice.Process')
    def test_start(self, mock_process, mock_connection):
        dice = Dice(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1],
            }
        )
        dice.start(1)

        self.assertTrue(mock_process.called)
        self.assertTrue(mock_connection.wait.called)

    def test_simulation_process_raises_exception(self):
        dice = Dice(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1],
            }
        )

        with self.assertRaises(SystemExit):
            dice._simulation_process(1, 1)

    @mock.patch('diceplayer.shared.external.dice.Dice._make_proc_dir')
    @mock.patch('diceplayer.shared.external.dice.Dice._make_dice_inputs')
    @mock.patch('diceplayer.shared.external.dice.Dice._run_dice')
    def test_simulation_process(self, mock_run_dice, mock_make_dice_inputs, mock_make_proc_dir):
        dice = Dice(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1],
            }
        )

        dice._simulation_process(1, 1)

        self.assertTrue(dice._make_proc_dir.called)
        self.assertTrue(dice._make_dice_inputs.called)
        self.assertTrue(dice._run_dice.called)


if __name__ == '__main__':
    unittest.main()
