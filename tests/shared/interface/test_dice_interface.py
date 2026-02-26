from diceplayer import logger
from diceplayer.config.player_config import PlayerConfig
from diceplayer.environment import Atom, Molecule, System
from diceplayer.interface import DiceInterface
from tests.mocks.mock_inputs import get_config_example
from tests.mocks.mock_proc import MockConnection, MockProc

import yaml

import io
import unittest
from unittest import mock


class TestDiceInterface(unittest.TestCase):
    def setUp(self):
        logger.set_logger(stream=io.StringIO())

        config = yaml.load(get_config_example(), Loader=yaml.Loader)
        self.config = PlayerConfig.model_validate(config["diceplayer"])

    def test_class_instantiation(self):
        dice = DiceInterface()

        self.assertIsInstance(dice, DiceInterface)

    def test_configure(self):
        dice = DiceInterface()

        self.assertIsNone(dice.step)
        self.assertIsNone(dice.system)

        # Ignoring the types for testing purposes
        dice.configure(self.config, System())

        self.assertIsNotNone(dice.step)
        self.assertIsNotNone(dice.system)

    def test_reset(self):
        dice = DiceInterface()

        dice.configure(self.config, System())

        self.assertTrue(hasattr(dice, "step"))
        self.assertTrue(hasattr(dice, "system"))

        dice.reset()

        self.assertFalse(hasattr(dice, "step"))
        self.assertFalse(hasattr(dice, "system"))

    @mock.patch("diceplayer.interface.dice_interface.Process", MockProc())
    @mock.patch("diceplayer.interface.dice_interface.connection", MockConnection)
    def test_start(self):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.start(1)

    @mock.patch("diceplayer.interface.dice_interface.connection", MockConnection)
    @mock.patch("diceplayer.interface.dice_interface.Process", MockProc(exitcode=1))
    def test_start_with_process_error(self):
        dice = DiceInterface()
        dice.configure(self.config, System())

        with self.assertRaises(SystemExit):
            dice.start(1)

    def test_simulation_process_raises_exception(self):
        dice = DiceInterface()

        with self.assertRaises(SystemExit):
            dice._simulation_process(1, 1)

    @mock.patch("diceplayer.interface.dice_interface.DiceInterface._make_proc_dir")
    @mock.patch("diceplayer.interface.dice_interface.DiceInterface._make_dice_inputs")
    @mock.patch("diceplayer.interface.dice_interface.DiceInterface._run_dice")
    def test_simulation_process(
        self, mock_run_dice, mock_make_dice_inputs, mock_make_proc_dir
    ):
        dice = DiceInterface()

        dice._simulation_process(1, 1)

        self.assertTrue(dice._make_proc_dir.called)
        self.assertTrue(dice._make_dice_inputs.called)
        self.assertTrue(dice._run_dice.called)

    @mock.patch("diceplayer.interface.dice_interface.Path.mkdir")
    @mock.patch("diceplayer.interface.dice_interface.Path.exists")
    def test_make_proc_dir_if_simdir_exists(self, mock_path_exists, mock_path_mkdir):
        dice = DiceInterface()
        dice.configure(self.config, System())

        mock_path_exists.return_value = False

        dice._make_proc_dir(1, 1)

        self.assertEqual(mock_path_mkdir.call_count, 2)

    @mock.patch("diceplayer.interface.dice_interface.Path.mkdir")
    @mock.patch("diceplayer.interface.dice_interface.Path.exists")
    def test_make_proc_dir_if_simdir_doesnt_exists(
        self, mock_path_exists, mock_path_mkdir
    ):
        dice = DiceInterface()
        dice.configure(self.config, System())

        mock_path_exists.return_value = False

        dice._make_proc_dir(1, 1)

        self.assertEqual(mock_path_mkdir.call_count, 2)

    def test_make_dice_seed(self):
        seed = DiceInterface._make_dice_seed()

        self.assertIsInstance(seed, int)

    def test_make_dice_inputs_nstep_len_two_with_randoninit_first_cycle_one(self):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.step.dice.nstep = [1, 1]

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

    @mock.patch("builtins.open", new_callable=mock.mock_open, read_data="test")
    @mock.patch("diceplayer.interface.dice_interface.Path.exists", return_value=True)
    def test_make_dice_inputs_nstep_len_two_with_randoninit_first_cycle_two(
        self, mock_path_exists, mock_open
    ):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.step.dice.nstep = [1, 1]

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

    @mock.patch("diceplayer.interface.dice_interface.Path.exists", return_value=False)
    def test_make_dice_inputs_raises_exception_on_last_not_found(
        self, mock_path_exists
    ):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.step.dice.nstep = [1, 1]

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
        dice = DiceInterface()
        dice.configure(self.config, System())

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

    @mock.patch("diceplayer.interface.dice_interface.os")
    @mock.patch("diceplayer.interface.dice_interface.shutil")
    @mock.patch("diceplayer.interface.dice_interface.Path.exists", return_value=True)
    def test_run_dice_on_first_cycle_run_successful(
        self, mock_path_exists, mock_shutils, mock_os
    ):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.step.dice.nstep = [1, 1, 1]

        dice.run_dice_file = mock.Mock()

        dice._run_dice(1, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 3)
        self.assertTrue(mock_shutils.copy.called)

        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.step.dice.nstep = [1, 1]

        dice.run_dice_file = mock.Mock()

        dice._run_dice(1, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 2)
        self.assertTrue(mock_shutils.copy.called)

    @mock.patch("diceplayer.interface.dice_interface.os")
    @mock.patch("diceplayer.interface.dice_interface.shutil")
    @mock.patch("diceplayer.interface.dice_interface.Path.exists", return_value=True)
    def test_run_dice_on_second_cycle_run_successful(
        self, mock_path_exists, mock_shutils, mock_os
    ):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.run_dice_file = mock.Mock()

        dice._run_dice(2, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 2)
        self.assertTrue(mock_shutils.copy.called)

        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.run_dice_file = mock.Mock()

        dice._run_dice(2, 1)

        self.assertTrue(mock_os.getcwd.called)
        self.assertTrue(mock_os.chdir.called)

        self.assertEqual(dice.run_dice_file.call_count, 2)
        self.assertTrue(mock_shutils.copy.called)

    @mock.patch("diceplayer.interface.dice_interface.os")
    @mock.patch("diceplayer.interface.dice_interface.shutil")
    @mock.patch("diceplayer.interface.dice_interface.Path.exists", return_value=False)
    def test_run_dice_raises_filenotfound_on_invalid_file(
        self, mock_path_exists, mock_shutils, mock_os
    ):
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.run_dice_file = mock.Mock()

        with self.assertRaises(FileNotFoundError):
            dice._run_dice(1, 1)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    def test_make_init_file(self, mock_open):
        example_atom = Atom(
            lbl=1,
            na=1,
            rx=1.0,
            ry=1.0,
            rz=1.0,
            chg=1.0,
            eps=1.0,
            sig=1.0,
        )

        main_molecule = Molecule("main_molecule")
        main_molecule.add_atom(example_atom)

        secondary_molecule = Molecule("secondary_molecule")
        secondary_molecule.add_atom(example_atom)

        system = System()
        system.add_type(main_molecule)
        system.add_type(secondary_molecule)

        dice = DiceInterface()
        dice.configure(self.config, system)

        dice.step.dice.nmol = [1, 1]

        last_xyz_file = io.StringIO()
        last_xyz_file.writelines(
            [
                "   TEST\n",
                " Configuration number :     TEST =   TEST  TEST  TEST\n",
                "   H     1.00000     1.00000     1.00000\n",
                "   H     1.00000     1.00000     1.00000\n",
            ]
        )
        last_xyz_file.seek(0)

        dice._make_init_file("test", last_xyz_file)

        mock_handler = mock_open()
        calls = mock_handler.write.call_args_list

        lines = list(map(lambda x: x[0][0], calls))

        expected_lines = [
            "  1.000000    1.000000    1.000000\n",
            "  1.000000    1.000000    1.000000\n",
            "$end",
        ]

        self.assertEqual(lines, expected_lines)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    def test_new_density(self, mock_open):
        example_atom = Atom(
            lbl=1,
            na=1,
            rx=1.0,
            ry=1.0,
            rz=1.0,
            chg=1.0,
            eps=1.0,
            sig=1.0,
        )

        main_molecule = Molecule("main_molecule")
        main_molecule.add_atom(example_atom)

        secondary_molecule = Molecule("secondary_molecule")
        secondary_molecule.add_atom(example_atom)

        system = System()
        system.add_type(main_molecule)
        system.add_type(secondary_molecule)

        dice = DiceInterface()
        dice.configure(self.config, system)

        last_xyz_file = io.StringIO()
        last_xyz_file.writelines(
            [
                "   TEST\n",
                " Configuration number :     TEST =   1  1  1\n",
                "   H     1.00000     1.00000     1.00000\n",
                "   H     1.00000     1.00000     1.00000\n",
            ]
        )
        last_xyz_file.seek(0)

        density = dice._new_density(last_xyz_file)

        self.assertEqual(density, 85.35451545000001)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    @mock.patch("diceplayer.interface.dice_interface.random")
    def test_make_nvt_ter(self, mock_random, mock_open):
        mock_random.random.return_value = 1

        dice = DiceInterface()
        dice.configure(self.config, System())

        dice._make_nvt_ter(1, "test")

        mock_handler = mock_open()
        calls = mock_handler.write.call_args_list

        lines = list(map(lambda x: x[0][0], calls))

        expected_lines = [
            "title = Diceplayer run - NVT Thermalization\n",
            "ncores = 4\n",
            "ljname = phb.ljc\n",
            "outname = phb\n",
            "nmol = 1 50\n",
            "dens = 0.75\n",
            "temp = 300.0\n",
            "init = yes\n",
            "nstep = 2000\n",
            "vstep = 0\n",
            "mstop = 1\n",
            "accum = no\n",
            "iprint = 1\n",
            "isave = 0\n",
            "irdf = 0\n",
            "seed = 1000000\n",
            "upbuf = 360",
        ]

        self.assertEqual(lines, expected_lines)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    @mock.patch("diceplayer.interface.dice_interface.random")
    def test_make_nvt_eq(self, mock_random, mock_open):
        mock_random.random.return_value = 1

        dice = DiceInterface()
        dice.configure(self.config, System())

        dice._make_nvt_eq(1, "test")

        mock_handler = mock_open()
        calls = mock_handler.write.call_args_list

        lines = list(map(lambda x: x[0][0], calls))

        expected_lines = [
            "title = Diceplayer run - NVT Production\n",
            "ncores = 4\n",
            "ljname = phb.ljc\n",
            "outname = phb\n",
            "nmol = 1 50\n",
            "dens = 0.75\n",
            "temp = 300.0\n",
            "init = no\n",
            "nstep = 3000\n",
            "vstep = 0\n",
            "mstop = 1\n",
            "accum = no\n",
            "iprint = 1\n",
            "isave = 1000\n",
            "irdf = 40\n",
            "seed = 1000000\n",
        ]

        self.assertEqual(lines, expected_lines)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    @mock.patch("diceplayer.interface.dice_interface.random")
    def test_make_npt_ter(self, mock_random, mock_open):
        mock_random.random.return_value = 1

        dice = DiceInterface()
        dice.configure(self.config, System())

        dice._make_npt_ter(1, "test")

        mock_handler = mock_open()
        calls = mock_handler.write.call_args_list

        lines = list(map(lambda x: x[0][0], calls))

        expected_lines = [
            "title = Diceplayer run - NPT Thermalization\n",
            "ncores = 4\n",
            "ljname = phb.ljc\n",
            "outname = phb\n",
            "nmol = 1 50\n",
            "press = 1.0\n",
            "temp = 300.0\n",
            "init = no\n",
            "vstep = 600\n",
            "nstep = 5\n",
            "mstop = 1\n",
            "accum = no\n",
            "iprint = 1\n",
            "isave = 0\n",
            "irdf = 0\n",
            "seed = 1000000\n",
        ]

        self.assertEqual(lines, expected_lines)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    @mock.patch("diceplayer.interface.dice_interface.random")
    def test_make_npt_eq(self, mock_random, mock_open):
        mock_random.random.return_value = 1

        dice = DiceInterface()
        dice.configure(self.config, System())

        dice._make_npt_eq("test")

        mock_handler = mock_open()
        calls = mock_handler.write.call_args_list

        lines = list(map(lambda x: x[0][0], calls))

        expected_lines = [
            "title = Diceplayer run - NPT Production\n",
            "ncores = 4\n",
            "ljname = phb.ljc\n",
            "outname = phb\n",
            "nmol = 1 50\n",
            "press = 1.0\n",
            "temp = 300.0\n",
            "nstep = 5\n",
            "vstep = 800\n",
            "init = no\n",
            "mstop = 1\n",
            "accum = no\n",
            "iprint = 1\n",
            "isave = 1000\n",
            "irdf = 40\n",
            "seed = 1000000\n",
        ]

        self.assertEqual(lines, expected_lines)

    @mock.patch("builtins.open", new_callable=mock.mock_open)
    def test_make_potentials(self, mock_open):
        example_atom = Atom(
            lbl=1,
            na=1,
            rx=1.0,
            ry=1.0,
            rz=1.0,
            chg=1.0,
            eps=1.0,
            sig=1.0,
        )

        main_molecule = Molecule("main_molecule")
        main_molecule.add_atom(example_atom)

        secondary_molecule = Molecule("secondary_molecule")
        secondary_molecule.add_atom(example_atom)

        system = System()
        system.add_type(main_molecule)
        system.add_type(secondary_molecule)

        dice = DiceInterface()
        dice.configure(self.config, system)

        dice._make_potentials("test")

        mock_handler = mock_open()
        calls = mock_handler.write.call_args_list

        lines = list(map(lambda x: x[0][0], calls))

        expected_lines = [
            "*\n",
            "2\n",
            "1 main_molecule\n",
            "1     1     1.00000    1.00000    1.00000    1.000000   1.00000  1.0000\n",
            "1 secondary_molecule\n",
            "1     1     1.00000    1.00000    1.00000    1.000000   1.00000  1.0000\n",
        ]

        self.assertEqual(lines, expected_lines)

    @mock.patch("diceplayer.interface.dice_interface.subprocess")
    @mock.patch(
        "builtins.open",
        new_callable=mock.mock_open,
        read_data="End of simulation\nBLABLA",
    )
    def test_run_dice_file(self, mock_open, mock_subprocess):
        mock_subprocess.call.return_value = 0
        dice = DiceInterface()
        dice.configure(self.config, System())

        dice.run_dice_file(1, 1, "test")

        self.assertTrue(mock_subprocess.call.called)
        self.assertTrue(mock_open.called)

    @mock.patch("diceplayer.interface.dice_interface.subprocess")
    @mock.patch("builtins.open", new_callable=mock.mock_open, read_data="Error\nBLABLA")
    def test_run_dice_file_raises_runtime_error_on_dice_file(
        self, mock_open, mock_subprocess
    ):
        mock_subprocess.call.return_value = 0
        dice = DiceInterface()
        dice.configure(self.config, System())

        with self.assertRaises(RuntimeError):
            dice.run_dice_file(1, 1, "test")

    @mock.patch("diceplayer.interface.dice_interface.subprocess")
    @mock.patch(
        "builtins.open",
        new_callable=mock.mock_open,
        read_data="End of simulation\nBLABLA",
    )
    def test_run_dice_file_raises_runtime_error_of_dice_exit_code(
        self, mock_open, mock_subprocess
    ):
        mock_subprocess.call.return_value = 1
        dice = DiceInterface()
        dice.configure(self.config, System())

        with self.assertRaises(RuntimeError):
            dice.run_dice_file(1, 1, "test")


if __name__ == "__main__":
    unittest.main()
