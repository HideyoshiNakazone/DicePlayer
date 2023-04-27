from diceplayer.shared.utils.logger import Logger

import unittest


class TestLogger(unittest.TestCase):
    def test_class_instantiation(self):
        logger = Logger()

        self.assertIsInstance(logger, Logger)

    def test_singleton(self):
        logger1 = Logger()
        logger2 = Logger()

        self.assertIs(logger1, logger2)

    def test_set_logger(self):
        logger = Logger()
        logger.set_logger('test_logger')

        self.assertIsNotNone(logger._logger)
        self.assertEqual(logger._logger.name, 'test_logger')

    def test_close(self):
        logger = Logger()
        logger.set_logger('test_logger')
        logger.close()

        self.assertEqual(len(logger._logger.handlers), 0)


if __name__ == '__main__':
    unittest.main()
