from diceplayer.environment import Atom, Molecule

import numpy as np
import numpy.testing as npt

import unittest


class TestMolecule(unittest.TestCase):
    def test_class_instantiation(self):
        mol = Molecule("test")

        self.assertIsInstance(mol, Molecule)

    def test_add_atom(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        self.assertEqual(len(mol.atom), 1)
        npt.assert_equal(mol.com, [1.0, 1.0, 1.0])

    def test_center_of_mass(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )
        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        npt.assert_equal(mol.com, [0.5, 0.5, 0.5])

    def test_center_of_mass_to_origin(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        mol.move_center_of_mass_to_origin()

        npt.assert_equal(mol.com, [0, 0, 0])

    def test_charges_and_dipole(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        actual_charge_dipole_array = mol.charges_and_dipole()

        expected_charge_dipole_array = [1.0, 0.0, 0.0, 0.0, 0.0]

        npt.assert_equal(actual_charge_dipole_array, expected_charge_dipole_array)

    def test_distances_between_atoms(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )
        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected = [
            [0.0, 1.73205081],
            [1.73205081, 0.0]
        ]
        actual = mol.distances_between_atoms()

        npt.assert_almost_equal(
            expected, actual
        )

    def test_inertia_tensor(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )
        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_inertia_tensor = [
            [1.00790, -0.50395, -0.50395],
            [-0.50395, 1.0079, -0.50395],
            [-0.50395, -0.50395, 1.0079],
        ]

        actual_inertia_tensor = mol.inertia_tensor

        npt.assert_equal(expected_inertia_tensor, actual_inertia_tensor)

    def test_principal_axes(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_evals, expected_evecs = (
            [0.0, 0.0, 0.0],
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        )

        evals, evecs = mol.principal_axes()

        npt.assert_equal(expected_evals, evals)
        npt.assert_equal(expected_evecs, evecs)

    def test_read_position(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_position = mol.read_position()

        actual_position = mol.read_position()

        npt.assert_equal(expected_position, actual_position)

    def test_update_charges(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_charges = [2.0]
        mol.update_charges(expected_charges)

        actual_charges = list(map(lambda a: a.chg, mol.atom))

        npt.assert_equal(expected_charges, actual_charges)

    def test_sizes_of_molecule(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        sizes = mol.sizes_of_molecule()

        expected_sizes = [0.0, 0.0, 0.0]

        npt.assert_equal(sizes, expected_sizes)

    def test_standard_orientation(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        mol.rotate_to_standard_orientation()

        expected_position = [0.0, 0.0, 0.0]

        self.assertEqual(mol.read_position().tolist(), expected_position)

    def test_translate(self):
        mol = Molecule("test")

        mol.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=1.0, rz=1.0, chg=1.0, eps=1.0, sig=1.0)
        )

        new_mol = mol.translate(np.array([-1, -1, -1]))

        expected_position = [0.0, 0.0, 0.0]

        self.assertEqual(new_mol.read_position().tolist(), expected_position)

    def test_minimum_distance(self):
        mol1 = Molecule("test1")
        mol1.add_atom(
            Atom(lbl=1, na=1, rx=0.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        mol2 = Molecule("test2")
        mol2.add_atom(
            Atom(lbl=1, na=1, rx=1.0, ry=0.0, rz=0.0, chg=1.0, eps=1.0, sig=1.0)
        )

        expected_distance = 1.0

        actual_distance = mol1.minimum_distance(mol2)

        self.assertEqual(expected_distance, actual_distance)


if __name__ == "__main__":
    unittest.main()
