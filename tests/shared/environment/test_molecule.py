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
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        self.assertEqual(len(mol.atom), 1)
        npt.assert_equal(mol.com, [1., 1., 1.])

    def test_center_of_mass(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )
        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        npt.assert_equal(mol.com, [.5, .5, .5])

    def test_center_of_mass_to_origin(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        mol.center_of_mass_to_origin()

        npt.assert_equal(mol.com, [0, 0, 0])

    def test_charges_and_dipole(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        actual_charge_dipole_array = mol.charges_and_dipole()

        expected_charge_dipole_array = [1.0, 0.0, 0.0, 0.0, 0.0]

        npt.assert_equal(
            actual_charge_dipole_array,
            expected_charge_dipole_array
        )

    def test_distances_between_atoms(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )
        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_distance_between_atoms = [[1.73205081], [1.73205081]]
        actual_distance_between_atoms = mol.distances_between_atoms()

        npt.assert_almost_equal(
            expected_distance_between_atoms,
            actual_distance_between_atoms
        )

    def test_inertia_tensor(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )
        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_inertia_tensor = [[1.00790, -0.50395, -0.50395],
                                   [-0.50395, 1.0079, -0.50395],
                                   [-0.50395, -0.50395, 1.0079]]

        actual_inertia_tensor = mol.inertia_tensor()

        npt.assert_equal(
            expected_inertia_tensor,
            actual_inertia_tensor
        )

    def test_principal_axes(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_evals, expected_evecs = [0., 0., 0.], [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]

        evals, evecs = mol.principal_axes()

        npt.assert_equal(expected_evals, evals)
        npt.assert_equal(expected_evecs, evecs)

    def test_read_position(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_position = mol.read_position()

        actual_position = mol.read_position()

        npt.assert_equal(
            expected_position,
            actual_position
        )

    def test_update_charges(self):
        mol = Molecule('test')

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_charges = [2.]
        mol.update_charges(expected_charges)

        actual_charges = list(map(lambda a: a.chg, mol.atom))

        npt.assert_equal(
            expected_charges,
            actual_charges
        )