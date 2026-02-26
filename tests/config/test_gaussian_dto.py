from diceplayer.config.gaussian_config import GaussianConfig

import unittest


class TestGaussianDTO(unittest.TestCase):
    def test_class_instantiation(self):
        gaussian_dto = GaussianConfig(
            level="test",
            qmprog="g16",
            keywords="test",
        )

        self.assertIsInstance(gaussian_dto, GaussianConfig)

    def test_is_valid_qmprog(self):
        with self.assertRaises(ValueError):
            gaussian_dto = GaussianConfig(
                level="test",
                qmprog="test",
                keywords="test",
            )

    def test_is_valid_level(self):
        with self.assertRaises(ValueError):
            gaussian_dto = GaussianConfig(
                level=None,
                qmprog="g16",
                keywords="test",
            )

    def test_from_dict(self):
        gaussian_dto = GaussianConfig.model_validate(
            {
                "level": "test",
                "qmprog": "g16",
                "keywords": "test",
            }
        )

        self.assertIsInstance(gaussian_dto, GaussianConfig)


if __name__ == "__main__":
    unittest.main()
