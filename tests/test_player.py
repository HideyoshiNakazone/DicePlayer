from diceplayer import VERSION, logger
from diceplayer.player import Player
from tests.mocks.mock_inputs import mock_open

import io
import unittest
from unittest import mock


class TestPlayer(unittest.TestCase):
    def setUp(self):
        logger.set_logger(stream=io.StringIO())

    @mock.patch("builtins.open", mock_open)
    def test_class_instantiation(self):
        # This file does not exist and it will be mocked
        player = Player("control.test.yml")

        self.assertIsInstance(player, Player)

    @mock.patch("builtins.open", mock_open)
    def test_start(self):
        player = Player("control.test.yml")

        player.gaussian_start = mock.MagicMock()
        player.dice_start = mock.MagicMock()

        player.start()

        self.assertEqual(player.dice_start.call_count, 3)
        self.assertEqual(player.gaussian_start.call_count, 3)

    @mock.patch("builtins.open", mock_open)
    @mock.patch("diceplayer.player.Path")
    def test_create_simulation_dir_if_already_exists(self, mock_path):
        player = Player("control.test.yml")
        mock_path.return_value.exists.return_value = True

        with self.assertRaises(FileExistsError):
            player.create_simulation_dir()

        self.assertTrue(mock_path.called)

    @mock.patch("builtins.open", mock_open)
    @mock.patch("diceplayer.player.Path")
    def test_create_simulation_dir_if_not_exists(self, mock_path):
        player = Player("control.test.yml")
        mock_path.return_value.exists.return_value = False

        player.create_simulation_dir()

        self.assertTrue(mock_path.called)

    @mock.patch("builtins.open", mock_open)
    @mock.patch("diceplayer.player.VERSION", "test")
    @mock.patch("diceplayer.player.sys")
    @mock.patch("diceplayer.player.weekday_date_time")
    def test_print_keywords(self, mock_date_func, mock_sys):
        player = Player("control.test.yml")

        mock_sys.version = "TEST"
        mock_date_func.return_value = "00 Test 0000 at 00:00:00"

        with self.assertLogs() as cm:
            player.print_keywords()

        expected_output = [
            "INFO:diceplayer:##########################################################################################\n#############               Welcome to DICEPLAYER version test                #############\n##########################################################################################\n",
            "INFO:diceplayer:Your python version is TEST\n",
            "INFO:diceplayer:Program started on 00 Test 0000 at 00:00:00\n",
            "INFO:diceplayer:Environment variables:",
            "INFO:diceplayer:OMP_STACKSIZE = Not set\n",
            "INFO:diceplayer:------------------------------------------------------------------------------------------\n                         DICE variables being used in this run:\n------------------------------------------------------------------------------------------\n",
            "INFO:diceplayer:combrule = *",
            "INFO:diceplayer:dens = 0.75",
            "INFO:diceplayer:isave = 1000",
            "INFO:diceplayer:ljname = phb.ljc",
            "INFO:diceplayer:nmol = [ 1 50 ]",
            "INFO:diceplayer:nstep = [ 2000 3000 4000 ]",
            "INFO:diceplayer:outname = phb",
            "INFO:diceplayer:press = 1.0",
            "INFO:diceplayer:progname = ~/.local/bin/dice",
            "INFO:diceplayer:randominit = first",
            "INFO:diceplayer:temp = 300.0",
            "INFO:diceplayer:upbuf = 360",
            "INFO:diceplayer:------------------------------------------------------------------------------------------\n                         GAUSSIAN variables being used in this run:\n------------------------------------------------------------------------------------------\n",
            "INFO:diceplayer:chg_tol = 0.01",
            "INFO:diceplayer:chgmult = [ 0 1 ]",
            "INFO:diceplayer:keywords = freq",
            "INFO:diceplayer:level = MP2/aug-cc-pVDZ",
            "INFO:diceplayer:pop = chelpg",
            "INFO:diceplayer:qmprog = g16",
            "INFO:diceplayer:\n",
        ]

        self.assertEqual(cm.output, expected_output)

    def test_validate_atom_dict(self):
        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 0,
                    "na": 1,
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid number of fields for site 1 for molecule type 1.",
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": "",
                    "na": 1,
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": 1.0,
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception), "Invalid lbl fields for site 1 for molecule type 1."
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": "",
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": 1.0,
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception), "Invalid na fields for site 1 for molecule type 1."
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": 1,
                    "rx": "",
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": 1.0,
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid rx fields for site 1 for molecule type 1. Value must be a float.",
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": 1,
                    "rx": 1.0,
                    "ry": "",
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": 1.0,
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid ry fields for site 1 for molecule type 1. Value must be a float.",
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": 1,
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": "",
                    "chg": 1.0,
                    "eps": 1.0,
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid rz fields for site 1 for molecule type 1. Value must be a float.",
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": 1,
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": "",
                    "eps": 1.0,
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid chg fields for site 1 for molecule type 1. Value must be a float.",
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": 1,
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": "",
                    "sig": 1.0,
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid eps fields for site 1 for molecule type 1. Value must be a float.",
        )

        with self.assertRaises(ValueError) as context:
            Player.validate_atom_dict(
                molecule_type=0,
                molecule_site=0,
                atom_dict={
                    "lbl": 1.0,
                    "na": 1,
                    "rx": 1.0,
                    "ry": 1.0,
                    "rz": 1.0,
                    "chg": 1.0,
                    "eps": 1.0,
                    "sig": "",
                },
            )
        self.assertEqual(
            str(context.exception),
            "Invalid sig fields for site 1 for molecule type 1. Value must be a float.",
        )

    @mock.patch("builtins.open", mock_open)
    @mock.patch("diceplayer.player.Path.exists", return_value=True)
    def test_read_potentials(self, mock_path_exists):
        player = Player("control.test.yml")

        player.read_potentials()

        self.assertEqual(player.system.molecule[0].molname, "TEST")
        self.assertEqual(len(player.system.molecule[0].atom), 1)

        self.assertEqual(player.system.molecule[1].molname, "PLACEHOLDER")
        self.assertEqual(len(player.system.molecule[1].atom), 1)

    @mock.patch("builtins.open", mock_open)
    @mock.patch("diceplayer.player.Path.exists")
    def test_read_potentials_error(self, mock_path_exists):
        player = Player("control.test.yml")

        # Testing file not found error
        mock_path_exists.return_value = False
        with self.assertRaises(RuntimeError) as context:
            player.read_potentials()

        self.assertEqual(str(context.exception), "Potential file phb.ljc not found.")

        # Enabling file found for next tests
        mock_path_exists.return_value = True

        # Testing combrule error
        with self.assertRaises(SystemExit) as context:
            player.config.dice.ljname = "phb.error.combrule.ljc"
            player.read_potentials()

        self.assertEqual(
            str(context.exception),
            "Error: expected a '*' or a '+' sign in 1st line of file phb.error.combrule.ljc",
        )

        # Testing ntypes error
        with self.assertRaises(SystemExit) as context:
            player.config.dice.ljname = "phb.error.ntypes.ljc"
            player.read_potentials()

        self.assertEqual(
            str(context.exception),
            "Error: expected an integer in the 2nd line of file phb.error.ntypes.ljc",
        )

        # Testing ntypes error on config
        with self.assertRaises(SystemExit) as context:
            player.config.dice.ljname = "phb.error.ntypes.config.ljc"
            player.read_potentials()

        self.assertEqual(
            str(context.exception),
            "Error: number of molecule types in file phb.error.ntypes.config.ljc "
            "must match that of 'nmol' keyword in config file",
        )

        # Testing nsite error
        with self.assertRaises(ValueError) as context:
            player.config.dice.ljname = "phb.error.nsites.ljc"
            player.read_potentials()

        self.assertEqual(
            str(context.exception),
            "Error: expected nsites to be an integer for molecule type 1",
        )

        # Testing molname error
        with self.assertRaises(ValueError) as context:
            player.config.dice.ljname = "phb.error.molname.ljc"
            player.read_potentials()

        self.assertEqual(
            str(context.exception),
            "Error: expected nsites and molname for the molecule type 1",
        )

    @mock.patch("builtins.open", mock_open)
    @mock.patch("diceplayer.player.Path.exists", return_value=True)
    def test_print_potentials(self, mock_path_exists):
        player = Player("control.test.yml")
        player.read_potentials()

        with self.assertLogs(level="INFO") as context:
            player.print_potentials()

        expected_output = [
            "INFO:diceplayer:==========================================================================================\n                    Potential parameters from file phb.ljc:\n------------------------------------------------------------------------------------------\n",
            "INFO:diceplayer:Combination rule: *",
            "INFO:diceplayer:Types of molecules: 2\n",
            "INFO:diceplayer:1 atoms in molecule type 1:",
            "INFO:diceplayer:---------------------------------------------------------------------------------",
            "INFO:diceplayer:Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass",
            "INFO:diceplayer:---------------------------------------------------------------------------------",
            "INFO:diceplayer:1     1     0.00000    0.00000    0.00000    0.000000   0.00000  0.0000    1.0079",
            "INFO:diceplayer:\n",
            "INFO:diceplayer:1 atoms in molecule type 2:",
            "INFO:diceplayer:---------------------------------------------------------------------------------",
            "INFO:diceplayer:Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass",
            "INFO:diceplayer:---------------------------------------------------------------------------------",
            "INFO:diceplayer:1     1     0.00000    0.00000    0.00000    0.000000   0.00000  0.0000    1.0079",
            "INFO:diceplayer:\n",
        ]

        self.assertEqual(context.output, expected_output)

    @mock.patch("builtins.open", mock_open)
    def test_dice_start(self):
        player = Player("control.test.yml")
        player.dice_interface = mock.MagicMock()
        player.dice_interface.start = mock.MagicMock()

        player.dice_start(1)

        player.dice_interface.start.assert_called_once()


if __name__ == "__main__":
    unittest.main()
