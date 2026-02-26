from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.environment.system import System

import unittest


class TestSystem(unittest.TestCase):
    def test_class_instantiation(self):
        system = System()

        self.assertIsInstance(system, System)

    def test_add_type(self):
        system = System()
        system.add_type(Molecule("test"))

        self.assertIsInstance(system.molecule, list)

        with self.assertRaises(TypeError) as ex:
            system.add_type("test")
            self.assertEqual(ex.exception, "Error: molecule is not a Molecule instance")


if __name__ == "__main__":
    unittest.main()
