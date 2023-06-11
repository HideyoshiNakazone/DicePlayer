from diceplayer.shared.interface.gaussian_interface import GaussianInterface
from diceplayer.shared.config.player_config import PlayerConfig
from diceplayer.shared.environment.system import System
from diceplayer import logger

from tests.mocks.mock_inputs import get_config_example

import yaml
import io


from unittest import mock
import unittest


class TestGaussianInterface(unittest.TestCase):
    def setUp(self) -> None:
        logger.set_logger(stream=io.StringIO())

        config = yaml.load(get_config_example(), Loader=yaml.Loader)
        self.config = PlayerConfig.from_dict(config['diceplayer'])

    def test_class_instantiation(self):
        gaussian_interface = GaussianInterface()
        self.assertIsInstance(gaussian_interface, GaussianInterface)

    def test_configure(self):
        gaussian_interface = GaussianInterface()

        self.assertIsNone(gaussian_interface.step)
        self.assertIsNone(gaussian_interface.system)

        gaussian_interface.configure(self.config, System())

        self.assertIsNotNone(gaussian_interface.step)
        self.assertIsNotNone(gaussian_interface.system)

    def test_reset(self):
        gaussian_interface = GaussianInterface()

        gaussian_interface.configure(self.config, System())

        self.assertIsNotNone(gaussian_interface.step)
        self.assertIsNotNone(gaussian_interface.system)

        gaussian_interface.reset()

        self.assertFalse(hasattr(gaussian_interface, 'step'))
        self.assertFalse(hasattr(gaussian_interface, 'system'))

    @mock.patch('diceplayer.shared.interface.gaussian_interface.Path.mkdir')
    @mock.patch('diceplayer.shared.interface.gaussian_interface.Path.exists')
    def test_make_qm_dir(self, mock_exists, mock_mkdir):
        mock_exists.return_value = False

        gaussian_interface = GaussianInterface()
        gaussian_interface.configure(self.config, System())

        gaussian_interface._make_qm_dir(1)

        mock_exists.assert_called_once()
        mock_mkdir.assert_called_once()

    @mock.patch('diceplayer.shared.interface.gaussian_interface.shutil.copy')
    @mock.patch('diceplayer.shared.interface.gaussian_interface.Path.exists')
    def test_copy_chk_file_from_previous_step(self, mock_exists, mock_copy):
        gaussian_interface = GaussianInterface()
        gaussian_interface.configure(self.config, System())

        mock_exists.side_effect = [False, True]

        gaussian_interface._copy_chk_file_from_previous_step(2)

        self.assertTrue(mock_exists.called)
        self.assertTrue(mock_copy.called)

    @mock.patch('diceplayer.shared.interface.gaussian_interface.shutil.copy')
    @mock.patch('diceplayer.shared.interface.gaussian_interface.Path.exists')
    def test_copy_chk_file_from_previous_step_no_previous_step(self, mock_exists, mock_copy):
        gaussian_interface = GaussianInterface()
        gaussian_interface.configure(self.config, System())

        mock_exists.side_effect = [False, False]

        with self.assertRaises(FileNotFoundError):
            gaussian_interface._copy_chk_file_from_previous_step(2)

    @mock.patch('diceplayer.shared.interface.gaussian_interface.shutil.copy')
    @mock.patch('diceplayer.shared.interface.gaussian_interface.Path.exists')
    def test_copy_chk_file_from_previous_step_current_exists(self, mock_exists, mock_copy):
        gaussian_interface = GaussianInterface()
        gaussian_interface.configure(self.config, System())

        mock_exists.side_effect = [True, True]

        with self.assertRaises(FileExistsError):
            gaussian_interface._copy_chk_file_from_previous_step(2)

    # def test_start(self):
    #     gaussian_interface = GaussianInterface()
    #     gaussian_interface.configure(self.config, System())
    #
    #     gaussian_interface._make_qm_dir = mock.Mock()
    #     gaussian_interface._copy_chk_file_from_previous_step = mock.Mock()
    #
    #     gaussian_interface.start(2)
    #
    #     gaussian_interface._make_qm_dir.assert_called_once_with(2)
    #     gaussian_interface._copy_chk_file_from_previous_step.assert_called_once_with(2)


if __name__ == '__main__':
    unittest.main()
