from diceplayer.shared.utils.logger import Logger, valid_logger

import logging

from unittest import mock
import unittest


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
        logger = Logger('test')

        self.assertIsInstance(logger, Logger)

    @mock.patch('builtins.open', mock.mock_open())
    def test_set_logger(self):
        logger = Logger('test')
        logger.set_logger()

        self.assertIsNotNone(logger._logger)
        self.assertEqual(logger._logger.name, 'test')

    @mock.patch('builtins.open', mock.mock_open())
    def test_close(self):
        logger = Logger('test')
        logger.set_logger()
        logger.close()

        self.assertEqual(len(logger._logger.handlers), 0)

    @mock.patch('builtins.open', mock.mock_open())
    def test_info(self):
        logger = Logger('test')
        logger.set_logger()

        with self.assertLogs(level='INFO') as cm:
            logger.info('test')

        self.assertEqual(cm.output, ['INFO:test:test'])

    @mock.patch('builtins.open', mock.mock_open())
    def test_debug(self):
        logger = Logger('test')
        logger.set_logger(level=logging.DEBUG)

        with self.assertLogs(level='DEBUG') as cm:
            logger.debug('test')

        self.assertEqual(cm.output, ['DEBUG:test:test'])

    @mock.patch('builtins.open', mock.mock_open())
    def test_warning(self):
        logger = Logger('test')
        logger.set_logger()

        with self.assertLogs(level='WARNING') as cm:
            logger.warning('test')

        self.assertEqual(cm.output, ['WARNING:test:test'])

    @mock.patch('builtins.open', mock.mock_open())
    def test_error(self):
        logger = Logger('test')
        logger.set_logger()

        with self.assertLogs(level='ERROR') as cm:
            logger.error('test')

        self.assertEqual(cm.output, ['ERROR:test:test'])


if __name__ == '__main__':
    unittest.main()
