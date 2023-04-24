from diceplayer.shared.environment.atom import Atom
from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.environment.system import System

import unittest


class TestSystem(unittest.TestCase):
    def test_class_instantiation(self):
        system = System()

        self.assertIsInstance(system, System)

    def test_add_type(self):
        system = System()
        system.add_type(0, Molecule('test'))

        self.assertIsInstance(system.molecule, list)
        self.assertIsInstance(system.nmols, list)

        with self.assertRaises(TypeError) as ex:
            system.add_type(0, 'test')
            self.assertEqual(ex.exception, 'Error: molecule is not a Molecule instance')

        with self.assertRaises(TypeError) as ex:
            system.add_type('test', Molecule('test'))
            self.assertEqual(ex.exception, 'Error: nmols is not an integer')


if __name__ == '__main__':
    unittest.main()
