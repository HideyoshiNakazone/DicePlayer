from diceplayer.utils import Logger, valid_logger

import io
import logging
import unittest
from unittest import mock


class TestValidateLogger(unittest.TestCase):
    def test_validate_logger(self):
        class MockLogger:
            _was_set = True

            @valid_logger
            def test_func(self):
                pass

        MockLogger().test_func()

    def test_validate_logger_exception(self):
        class MockLogger:
            _was_set = False

            @valid_logger
            def test_func(self):
                pass

        with self.assertRaises(AssertionError):
            MockLogger().test_func()


class TestLogger(unittest.TestCase):
    def test_class_instantiation(self):
        logger = Logger("test")

        self.assertIsInstance(logger, Logger)

    @mock.patch("builtins.open", mock.mock_open())
    def test_set_logger_to_file(self):
        logger = Logger("test")

        logger.set_logger(stream=io.StringIO())

        self.assertIsNotNone(logger._logger)
        self.assertEqual(logger._logger.name, "test")

    def test_set_logger_to_stream(self):
        logger = Logger("test")

        logger.set_logger(stream=io.StringIO())

        self.assertIsNotNone(logger._logger)
        self.assertEqual(logger._logger.name, "test")

    @mock.patch("builtins.open", mock.mock_open())
    @mock.patch("diceplayer.utils.logger.Path.exists")
    @mock.patch("diceplayer.utils.logger.Path.rename")
    def test_set_logger_if_file_exists(self, mock_rename, mock_exists):
        logger = Logger("test")

        mock_exists.return_value = True
        logger.set_logger()

        self.assertTrue(mock_rename.called)
        self.assertIsNotNone(logger._logger)
        self.assertEqual(logger._logger.name, "test")

    @mock.patch("builtins.open", mock.mock_open())
    @mock.patch("diceplayer.utils.logger.Path.exists")
    @mock.patch("diceplayer.utils.logger.Path.rename")
    def test_set_logger_if_file_not_exists(self, mock_rename, mock_exists):
        logger = Logger("test")

        mock_exists.return_value = False
        logger.set_logger()

        self.assertFalse(mock_rename.called)
        self.assertIsNotNone(logger._logger)
        self.assertEqual(logger._logger.name, "test")

    @mock.patch("builtins.open", mock.mock_open())
    def test_close(self):
        logger = Logger("test")

        logger.set_logger()
        logger.close()

        self.assertEqual(len(logger._logger.handlers), 0)

    @mock.patch("builtins.open", mock.mock_open())
    def test_info(self):
        logger = Logger("test")
        logger.set_logger(stream=io.StringIO())

        with self.assertLogs(level="INFO") as cm:
            logger.info("test")

        self.assertEqual(cm.output, ["INFO:test:test"])

    @mock.patch("builtins.open", mock.mock_open())
    def test_debug(self):
        logger = Logger("test")
        logger.set_logger(stream=io.StringIO(), level=logging.DEBUG)

        with self.assertLogs(level="DEBUG") as cm:
            logger.debug("test")

        self.assertEqual(cm.output, ["DEBUG:test:test"])

    @mock.patch("builtins.open", mock.mock_open())
    def test_warning(self):
        logger = Logger("test")
        logger.set_logger(stream=io.StringIO())

        with self.assertLogs(level="WARNING") as cm:
            logger.warning("test")

        self.assertEqual(cm.output, ["WARNING:test:test"])

    @mock.patch("builtins.open", mock.mock_open())
    def test_error(self):
        logger = Logger("test")
        logger.set_logger(stream=io.StringIO())

        with self.assertLogs(level="ERROR") as cm:
            logger.error("test")

        self.assertEqual(cm.output, ["ERROR:test:test"])


if __name__ == "__main__":
    unittest.main()
