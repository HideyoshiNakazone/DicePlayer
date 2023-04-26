from diceplayer.shared.interface.dice_interface import DiceInterface
from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.config.step_dto import StepDTO

import io

from tests.mocks.mock_proc import MockConnection, MockProc

from unittest import mock
import unittest


class TestDiceInterface(unittest.TestCase):
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

        self.assertIsNone(dice.step)

        dice.configure('test')

        self.assertIsNotNone(dice.step)

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

    @mock.patch('diceplayer.shared.interface.dice_interface.Process', MockProc())
    @mock.patch('diceplayer.shared.interface.dice_interface.connection', MockConnection)
    def test_start(self):
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

    @mock.patch('diceplayer.shared.interface.dice_interface.connection', MockConnection)
    @mock.patch('diceplayer.shared.interface.dice_interface.Process', MockProc(exitcode=1))
    def test_start_with_process_error(self):
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
                nprocs=2,
                simulation_dir='test',
                altsteps=1,
                molecule=[],
                nmol=[],
            )
        )

        with self.assertRaises(SystemExit):
            dice.start(1)

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

    @mock.patch('diceplayer.shared.interface.dice_interface.DiceInterface._make_proc_dir')
    @mock.patch('diceplayer.shared.interface.dice_interface.DiceInterface._make_dice_inputs')
    @mock.patch('diceplayer.shared.interface.dice_interface.DiceInterface._run_dice')
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

    @mock.patch('diceplayer.shared.interface.dice_interface.Path.mkdir')
    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists')
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

    @mock.patch('diceplayer.shared.interface.dice_interface.Path.mkdir')
    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists')
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

    def test_make_dice_inputs_nstep_len_two_with_randoninit_first_cycle_one(self):
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

        dice._make_potentials = mock.Mock()

        dice._make_init_file = mock.Mock()
        dice._new_density = mock.Mock()

        dice._make_nvt_ter = mock.Mock()
        dice._make_nvt_eq = mock.Mock()
        dice._make_npt_ter = mock.Mock()
        dice._make_npt_eq = mock.Mock()

        dice._make_dice_inputs(1, 1)

        self.assertTrue(dice._make_potentials.called)

        self.assertFalse(dice._make_init_file.called)
        self.assertFalse(dice._new_density.called)

        self.assertTrue(dice._make_nvt_ter.called)
        self.assertTrue(dice._make_nvt_eq.called)

        self.assertFalse(dice._make_npt_ter.called)
        self.assertFalse(dice._make_npt_eq.called)

    @mock.patch('builtins.open', new_callable=mock.mock_open, read_data='test')
    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists', return_value=True)
    def test_make_dice_inputs_nstep_len_two_with_randoninit_first_cycle_two(self, mock_path_exists, mock_open):
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

        dice._make_potentials = mock.Mock()

        dice._make_init_file = mock.Mock()
        dice._new_density = mock.Mock()

        dice._make_nvt_ter = mock.Mock()
        dice._make_nvt_eq = mock.Mock()
        dice._make_npt_ter = mock.Mock()
        dice._make_npt_eq = mock.Mock()

        dice._make_dice_inputs(2, 1)

        self.assertTrue(dice._make_potentials.called)

        self.assertTrue(dice._make_init_file.called)
        self.assertTrue(dice._new_density.called)

        self.assertFalse(dice._make_nvt_ter.called)
        self.assertTrue(dice._make_nvt_eq.called)

        self.assertFalse(dice._make_npt_ter.called)
        self.assertFalse(dice._make_npt_eq.called)

    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists', return_value=False)
    def test_make_dice_inputs_raises_exception_on_last_not_found(self, mock_path_exists):
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

        dice._make_potentials = mock.Mock()

        dice._make_init_file = mock.Mock()
        dice._new_density = mock.Mock()

        dice._make_nvt_ter = mock.Mock()
        dice._make_nvt_eq = mock.Mock()
        dice._make_npt_ter = mock.Mock()
        dice._make_npt_eq = mock.Mock()

        with self.assertRaises(FileNotFoundError):
            dice._make_dice_inputs(2, 1)

    def test_make_dice_inputs_nstep_len_three_with_randoninit_first_cycle_one(self):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1, 1],
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

        dice._make_potentials = mock.Mock()

        dice._make_init_file = mock.Mock()
        dice._new_density = mock.Mock()

        dice._make_nvt_ter = mock.Mock()
        dice._make_nvt_eq = mock.Mock()
        dice._make_npt_ter = mock.Mock()
        dice._make_npt_eq = mock.Mock()

        dice._make_dice_inputs(1, 1)

        self.assertTrue(dice._make_potentials.called)

        self.assertFalse(dice._make_init_file.called)
        self.assertFalse(dice._new_density.called)

        self.assertTrue(dice._make_nvt_ter.called)
        self.assertFalse(dice._make_nvt_eq.called)

        self.assertTrue(dice._make_npt_ter.called)
        self.assertTrue(dice._make_npt_eq.called)

    @mock.patch('diceplayer.shared.interface.dice_interface.os')
    @mock.patch('diceplayer.shared.interface.dice_interface.shutil')
    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists', return_value=True)
    def test_run_dice_on_first_cycle_run_successful(self, mock_path_exists, mock_shutils, mock_os):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1, 1],
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

        dice.run_dice_file = mock.Mock()

        dice._run_dice(1, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 3)
        self.assertTrue(mock_shutils.copy.called)

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

        dice.run_dice_file = mock.Mock()

        dice._run_dice(1, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 2)
        self.assertTrue(mock_shutils.copy.called)

    @mock.patch('diceplayer.shared.interface.dice_interface.os')
    @mock.patch('diceplayer.shared.interface.dice_interface.shutil')
    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists', return_value=True)
    def test_run_dice_on_second_cycle_run_successful(self, mock_path_exists, mock_shutils, mock_os):
        dice = DiceInterface(
            {
                'ljname': 'test',
                'outname': 'test',
                'ncores': 1,
                'dens': 1.0,
                'nmol': [1],
                'nstep': [1, 1, 1],
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

        dice.run_dice_file = mock.Mock()

        dice._run_dice(2, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 2)
        self.assertTrue(mock_shutils.copy.called)

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

        dice.run_dice_file = mock.Mock()

        dice._run_dice(2, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 1)
        self.assertTrue(mock_shutils.copy.called)

    @mock.patch('diceplayer.shared.interface.dice_interface.os')
    @mock.patch('diceplayer.shared.interface.dice_interface.shutil')
    @mock.patch('diceplayer.shared.interface.dice_interface.Path.exists', return_value=False)
    def test_run_dice_on_second_cycle_run_successful(self, mock_path_exists, mock_shutils, mock_os):
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

        dice.run_dice_file = mock.Mock()

        with self.assertRaises(FileNotFoundError):
            dice._run_dice(1, 1)

    def test_make_init_file(self):
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
                molecule=[
                    
                ],
                nmol=[],
            )
        )

        last_xyz_file = io.StringIO()
        last_xyz_file.writelines([
            'TEST',
            'TEST',
            'TEST',
            'TEST',
        ])


if __name__ == '__main__':
    unittest.main()
