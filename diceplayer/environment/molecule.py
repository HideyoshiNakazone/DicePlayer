from __future__ import annotations

from diceplayer import logger
from diceplayer.environment import Atom
from diceplayer.utils.cache import invalidate_computed_properties
from diceplayer.utils.misc import BOHR2ANG, EA_2_DEBYE
from diceplayer.utils.ptable import GHOST_NUMBER

import numpy as np
import numpy.typing as npt
from numpy.linalg import linalg
from typing_extensions import List, Self, Tuple

import math
from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property


@dataclass
class Molecule:
    """
    Molecule class declaration. This class is used throughout the DicePlayer program to represent molecules.

    Atributes:
        molname (str): The name of the represented molecule
        atom (List[Atom]): List of atoms of the represented molecule
        total_mass (int): The total mass of the molecule
        com (npt.NDArray[np.float64]): The center of mass of the molecule
        inertia_tensor (npt.NDArray[np.float64]): The inertia tensor of the molecule
    """

    molname: str
    atom: List[Atom] = field(default_factory=list)

    @cached_property
    def total_mass(self) -> float:
        return sum(atom.mass for atom in self.atom)

    @cached_property
    def com(self) -> npt.NDArray[np.float64]:
        com = np.zeros(3)

        for atom in self.atom:
            com += atom.mass * np.array([atom.rx, atom.ry, atom.rz])

        com = com / self.total_mass

        return com

    @cached_property
    def inertia_tensor(self) -> npt.NDArray[np.float64]:
        """
        Calculates the inertia tensor of the molecule.

        Returns:
            npt.NDArray[np.float64]: inertia tensor of the molecule.
        """
        inertia_tensor = np.zeros((3, 3), dtype=np.float64)

        for atom in self.atom:
            dx = atom.rx - self.com[0]
            dy = atom.ry - self.com[1]
            dz = atom.rz - self.com[2]

            inertia_tensor[0, 0] += atom.mass * (dy**2 + dz**2)
            inertia_tensor[1, 1] += atom.mass * (dz**2 + dx**2)
            inertia_tensor[2, 2] += atom.mass * (dx**2 + dy**2)

            inertia_tensor[0, 1] -= atom.mass * dx * dy
            inertia_tensor[0, 2] -= atom.mass * dx * dz
            inertia_tensor[1, 2] -= atom.mass * dy * dz

        # enforce symmetry
        inertia_tensor[1, 0] = inertia_tensor[0, 1]
        inertia_tensor[2, 0] = inertia_tensor[0, 2]
        inertia_tensor[2, 1] = inertia_tensor[1, 2]

        return inertia_tensor

    @invalidate_computed_properties()
    def add_atom(self, a: Atom) -> None:
        """
        Adds Atom instance to the molecule.

        Args:
            a (Atom): Atom instance to be added to atom list.
        """

        self.atom.append(a)

    @invalidate_computed_properties()
    def remove_atom(self, a: Atom) -> None:
        """
        Removes Atom instance from the molecule.

        Args:
            a (Atom): Atom instance to be removed from atom list.
        """

        self.atom.remove(a)

    @invalidate_computed_properties()
    def move_center_of_mass_to_origin(self) -> None:
        """
        Updated positions based on the center of mass of the molecule
        """
        for atom in self.atom:
            atom.rx -= self.com[0]
            atom.ry -= self.com[1]
            atom.rz -= self.com[2]

    @invalidate_computed_properties()
    def rotate_to_standard_orientation(self) -> None:
        """
        Rotates the molecule to the standard orientation
        """

        self.move_center_of_mass_to_origin()
        evals, evecs = self.principal_axes()

        if np.isclose(linalg.det(evecs), -1):
            evecs[:, 2] *= -1

        if not np.isclose(linalg.det(evecs), 1):
            raise RuntimeError(
                "Error: could not make a rotation matrix while adopting the standard orientation"
            )

        coords = np.array([(a.rx, a.ry, a.rz) for a in self.atom])
        rotated = coords @ evecs.T

        for atom, pos in zip(self.atom, rotated):
            atom.rx, atom.ry, atom.rz = pos

    def charges_and_dipole(self) -> List[float]:
        """
        Calculates the charges and dipole of the molecule atoms

        Returns:
            List[float]: Respectivly magnitude of the: charge magnitude, first dipole,
            second dipole, third dipole and total dipole.
        """

        charge = 0
        dipole = np.zeros(3)
        for atom in self.atom:
            position = np.array([atom.rx, atom.ry, atom.rz])
            dipole += atom.chg * position
            charge += atom.chg

        dipole *= EA_2_DEBYE
        total_dipole = math.sqrt(dipole[0] ** 2 + dipole[1] ** 2 + dipole[2] ** 2)

        return [charge, dipole[0], dipole[1], dipole[2], total_dipole]

    def distances_between_atoms(self) -> npt.NDArray[np.float64]:
        """
        Calculates distances between the atoms of the molecule

        Returns:
           NDArray[Shape["Any,Any"],Float]: distances between the atoms.
        """
        coords = np.array([(a.rx, a.ry, a.rz) for a in self.atom], dtype=np.float64)
        diff = coords[:, None, :] - coords[None, :, :]
        return np.linalg.norm(diff, axis=-1)

    def principal_axes(self) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """
        Calculates the principal axes of the molecule

        Returns:
            Tuple[np.ndarray, np.ndarray]: Tuple where the first value is the Eigen Values and the second is the Eigen Vectors,
             representing the principal axes of the molecule.
        """

        try:
            evals, evecs = linalg.eigh(self.inertia_tensor)

            idx = np.argsort(evals)
            evals = evals[idx]
            evecs = evecs[:, idx]
        except ValueError:
            raise RuntimeError(
                "Error: diagonalization of inertia tensor did not converge"
            )

        return evals, evecs

    def read_position(self) -> npt.NDArray[np.float64]:
        """Reads the position of the molecule from the position values of the atoms

        Returns:
            np.ndarray: internal position relative to atoms of the molecule
        """
        coords = np.array([(a.rx, a.ry, a.rz) for a in self.atom], dtype=np.float64)
        return coords.ravel() * BOHR2ANG

    def update_charges(self, charges: npt.NDArray[np.float64]) -> int:
        """
        Updates the charges of the atoms of the molecule and
        returns the max difference between the new and old charges
        """
        diff = 0
        for i, atom in enumerate(self.atom):
            diff = max(diff, abs(atom.chg - charges[i]))
            atom.chg = charges[i]

        return diff

    def sizes_of_molecule(self) -> List[float]:
        """
        Calculates sides of the smallest box that the molecule could fit

        Returns:
            List[float]: list of the sizes of the molecule
        """
        coords = np.array([(a.rx, a.ry, a.rz) for a in self.atom], dtype=np.float64)
        return (coords.max(axis=0) - coords.min(axis=0)).tolist()

    def translate(self, vector: np.ndarray) -> Self:
        """
        Creates a new Molecule object where its' atoms has been translated by a vector

        Args:
            vector (np.ndarray): translation vector

        Returns:
            Molecule: new Molecule object translated by a vector
        """
        vec = np.asarray(vector, dtype=np.float64)
        if vec.shape != (3,):
            raise ValueError("translation vector must be shape (3,)")

        new_molecule = deepcopy(self)

        for atom in new_molecule.atom:
            atom.rx += vector[0]
            atom.ry += vector[1]
            atom.rz += vector[2]

        return new_molecule

    def print_mol_info(self) -> None:
        """
        Prints the Molecule information into a Output File
        """

        logger.info(
            "    Center of mass = ( {:>10.4f} , {:>10.4f} , {:>10.4f} )".format(
                self.com[0], self.com[1], self.com[2]
            )
        )
        evals, evecs = self.principal_axes()

        logger.info(
            "    Moments of inertia =  {:>9E}  {:>9E}  {:>9E}".format(
                evals[0], evals[1], evals[2]
            )
        )

        logger.info(
            "    Major principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )".format(
                evecs[0, 0], evecs[1, 0], evecs[2, 0]
            )
        )
        logger.info(
            "    Inter principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )".format(
                evecs[0, 1], evecs[1, 1], evecs[2, 1]
            )
        )
        logger.info(
            "    Minor principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )".format(
                evecs[0, 2], evecs[1, 2], evecs[2, 2]
            )
        )

        sizes = self.sizes_of_molecule()
        logger.info(
            "    Characteristic lengths = ( {:>6.2f} , {:>6.2f} , {:>6.2f} )".format(
                sizes[0], sizes[1], sizes[2]
            )
        )
        logger.info("    Total mass = {:>8.2f} au".format(self.total_mass))

        chg_dip = self.charges_and_dipole()
        logger.info("    Total charge = {:>8.4f} e".format(chg_dip[0]))
        logger.info(
            "    Dipole moment = ( {:>9.4f} , {:>9.4f} , {:>9.4f} )     Total = {:>9.4f} Debye".format(
                chg_dip[1], chg_dip[2], chg_dip[3], chg_dip[4]
            )
        )

    def minimum_distance(self, molec: Self) -> float:
        """
        Return the minimum distance between two molecules

        Args:
            molec (Molecule): Molecule object to be compared

        Returns:
            float: minimum distance between the two molecules
        """
        coords_a = np.array(
            [(a.rx, a.ry, a.rz) for a in self.atom if a.na != GHOST_NUMBER]
        )
        coords_b = np.array(
            [(a.rx, a.ry, a.rz) for a in molec.atom if a.na != GHOST_NUMBER]
        )

        if len(coords_a) == 0 or len(coords_b) == 0:
            raise ValueError("No real atoms to compare")

        diff = coords_a[:, None, :] - coords_b[None, :, :]
        d2 = np.sum(diff**2, axis=-1)
        return np.sqrt(d2.min())
