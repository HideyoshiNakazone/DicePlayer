from diceplayer.shared.config.gaussian_config import GaussianDTO

import unittest


class TestGaussianDTO(unittest.TestCase):
    def test_class_instantiation(self):
        gaussian_dto = GaussianDTO(
            level='test',
            qmprog='g16',
            keywords='test',
        )

        self.assertIsInstance(gaussian_dto, GaussianDTO)

    def test_is_valid_qmprog(self):
        with self.assertRaises(ValueError):
            gaussian_dto = GaussianDTO(
                level='test',
                qmprog='test',
                keywords='test',
            )

    def test_is_valid_level(self):
        with self.assertRaises(ValueError):
            gaussian_dto = GaussianDTO(
                level=None,
                qmprog='g16',
                keywords='test',
            )

    def test_from_dict(self):
        gaussian_dto = GaussianDTO.from_dict(
            {
                'level': 'test',
                'qmprog': 'g16',
                'keywords': 'test',
            }
        )

        self.assertIsInstance(gaussian_dto, GaussianDTO)


if __name__ == '__main__':
    unittest.main()
