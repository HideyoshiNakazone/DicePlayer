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

    def test_center_of_mass_distance(self):
        system = System()

        a = Molecule('test')
        a.add_atom(
            Atom(lbl=0, na=1, rx=0, ry=0, rz=0, chg=0, eps=0, sig=0)
        )
        system.add_type(1, a)

        b = Molecule('test')
        b.add_atom(
            Atom(lbl=0, na=1, rx=0, ry=0, rz=0, chg=0, eps=0, sig=0)
        )
        system.add_type(1, b)

        self.assertIsInstance(system.center_of_mass_distance(0, 1), float)




if __name__ == '__main__':
    unittest.main()
