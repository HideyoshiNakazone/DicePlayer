from diceplayer.shared.config.step_dto import StepDTO
from diceplayer.shared.interface.dice_interface import DiceInterface

from unittest import mock
import unittest


class TestDice(unittest.TestCase):
    def test_class_instantiation(self):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )

        self.assertIsInstance(dice, DiceInterface)

    def test_configure(self):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )

        self.assertFalse(hasattr(dice, 'step'))

        dice.configure('test')

        self.assertTrue(hasattr(dice, 'step'))

    def test_reset(self):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )

        dice.configure('test')

        self.assertTrue(hasattr(dice, 'step'))

        dice.reset()

        self.assertFalse(hasattr(dice, 'step'))

    @mock.patch('diceplayer.shared.external.dice.connection')
    @mock.patch('diceplayer.shared.external.dice.Process')
    def test_start(self, mock_process, mock_connection):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )
        dice.configure(
            StepDTO(
                ncores=1,
                nprocs=1,
                simulation_dir='test',
                altsteps=1,
                molecule=[],
                nmol=[],
            )
        )

        dice.start(1)

        self.assertTrue(mock_process.called)
        self.assertTrue(mock_connection.wait.called)

    def test_simulation_process_raises_exception(self):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )

        with self.assertRaises(SystemExit):
            dice._simulation_process(1, 1)

    @mock.patch('diceplayer.shared.external.dice.Dice._make_proc_dir')
    @mock.patch('diceplayer.shared.external.dice.Dice._make_dice_inputs')
    @mock.patch('diceplayer.shared.external.dice.Dice._run_dice')
    def test_simulation_process(self, mock_run_dice, mock_make_dice_inputs, mock_make_proc_dir):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )

        dice._simulation_process(1, 1)

        self.assertTrue(dice._make_proc_dir.called)
        self.assertTrue(dice._make_dice_inputs.called)
        self.assertTrue(dice._run_dice.called)

    @mock.patch('diceplayer.shared.external.dice.Path.mkdir')
    @mock.patch('diceplayer.shared.external.dice.Path.exists')
    def test_make_proc_dir_if_simdir_exists(self, mock_path_exists, mock_path_mkdir):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )
        dice.configure(
            StepDTO(
                ncores=1,
                nprocs=1,
                simulation_dir='test',
                altsteps=1,
                molecule=[],
                nmol=[],
            )
        )

        mock_path_exists.return_value = False

        dice._make_proc_dir(1, 1)

        self.assertEqual(mock_path_mkdir.call_count, 2)

    @mock.patch('diceplayer.shared.external.dice.Path.mkdir')
    @mock.patch('diceplayer.shared.external.dice.Path.exists')
    def test_make_proc_dir_if_simdir_doesnt_exists(self, mock_path_exists, mock_path_mkdir):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )
        dice.configure(
            StepDTO(
                ncores=1,
                nprocs=1,
                simulation_dir='test',
                altsteps=1,
                molecule=[],
                nmol=[],
            )
        )

        mock_path_exists.return_value = False

        dice._make_proc_dir(1, 1)

        self.assertEqual(mock_path_mkdir.call_count, 2)

    def test_make_dice_seed(self):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1],
            }
        )

        seed = dice._make_dice_seed()

        self.assertIsInstance(seed, int)


if __name__ == '__main__':
    unittest.main()
