from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.environment.atom import Atom

import numpy.testing as npt
import unittest


class TestMolecule(unittest.TestCase):
    def test_class_instantiation(self):
        mol = Molecule('test')

        self.assertIsInstance(mol, Molecule)

    def test_add_atom(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(
                lbl=1,
                na=1,
                rx=1.0,
                ry=1.0,
                rz=1.0,
                chg=1.0,
                eps=1.0,
                sig=1.0,
            )
        )

        self.assertEqual(len(mol.atom), 1)
        npt.assert_equal(mol.com, [1., 1., 1.])

    def test_center_of_mass(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(
                lbl=1,
                na=1,
                rx=1.0,
                ry=1.0,
                rz=1.0,
                chg=1.0,
                eps=1.0,
                sig=1.0,
            )
        )
        mol.add_atom(
            Atom(
                lbl=1,
                na=1,
                rx=0.0,
                ry=0.0,
                rz=0.0,
                chg=1.0,
                eps=1.0,
                sig=1.0,
            )
        )
        
        npt.assert_equal(mol.com, [.5, .5, .5])
