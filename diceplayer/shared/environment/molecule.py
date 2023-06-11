from __future__ import annotations

from diceplayer import logger
from diceplayer.shared.environment.atom import Atom
from diceplayer.shared.utils.misc import BOHR2ANG
from diceplayer.shared.utils.ptable import ghost_number

import numpy as np
from nptyping import Float, NDArray, Shape
from numpy.linalg import linalg

import logging
import math
from copy import deepcopy
from typing import Any, List, Tuple, Union


class Molecule:
    """
    Molecule class declaration. This class is used throughout the DicePlayer program to represent molecules.

    Atributes:
        molname (str): The name of the represented molecule
        atom (List[Atom]): List of atoms of the represented molecule
        position (NDArray[Any, Any]): The position relative to the internal atoms of the represented molecule
        energy (NDArray[Any, Any]): The energy of the represented molecule
        gradient (NDArray[Any, Any]): The first derivative of the energy relative to the position
        hessian (NDArray[Any, Any]): The second derivative of the energy relative to the position
        total_mass (int): The total mass of the molecule
        com (NDArray[Any, Any]): The center of mass of the molecule
    """

    def __init__(self, molname: str) -> None:
        """
        The constructor function __init__ is used to create new instances of the Molecule class.

        Args:
            molname (str): Molecule name
        """
        self.molname: str = molname

        self.atom: List[Atom] = []
        self.position: NDArray[Any, Any]
        self.energy: NDArray[Any, Any]
        self.gradient: NDArray[Any, Any]
        self.hessian: NDArray[Any, Any]

        self.ghost_atoms: List[Atom] = []
        self.lp_atoms: List[Atom] = []

        self.total_mass: int = 0
        self.com: Union[None, NDArray[Any, Any]] = None

    def add_atom(self, a: Atom) -> None:
        """
        Adds Atom instance to the molecule.

        Args:
            a (Atom): Atom instance to be added to atom list.
        """

        self.atom.append(a)
        self.total_mass += a.mass

        self.center_of_mass()

    def center_of_mass(self) -> NDArray[Any, Any]:
        """
        Calculates the center of mass of the molecule
        """

        self.com = np.zeros(3)

        for atom in self.atom:
            self.com += atom.mass * np.array([atom.rx, atom.ry, atom.rz])

        self.com = self.com / self.total_mass

        return self.com

    def center_of_mass_to_origin(self) -> None:
        """
        Updated positions based on the center of mass of the molecule
        """
        for atom in self.atom:
            atom.rx -= self.com[0]
            atom.ry -= self.com[1]
            atom.rz -= self.com[2]

        self.center_of_mass()

    def charges_and_dipole(self) -> List[float]:
        """
        Calculates the charges and dipole of the molecule atoms

        Returns:
            List[float]: Respectivly magnitude of the: charge magnitude, first dipole,
            second dipole, third dipole and total dipole.
        """

        eA_to_Debye = 1 / 0.20819434
        charge = 0
        dipole = np.zeros(3)
        for atom in self.atom:
            position = np.array([atom.rx, atom.ry, atom.rz])
            dipole += atom.chg * position
            charge += atom.chg

        dipole *= eA_to_Debye
        total_dipole = math.sqrt(dipole[0] ** 2 + dipole[1] ** 2 + dipole[2] ** 2)

        return [charge, dipole[0], dipole[1], dipole[2], total_dipole]

    def distances_between_atoms(self) -> NDArray[Shape["Any,Any"], Float]:
        """
        Calculates distances between the atoms of the molecule

        Returns:
           NDArray[Shape["Any,Any"],Float]: distances between the atoms.
        """

        distances = []
        dim = len(self.atom)
        for index1, atom1 in enumerate(self.atom):
            for index2, atom2 in enumerate(self.atom):
                if index1 != index2:
                    dx = atom1.rx - atom2.rx
                    dy = atom1.ry - atom2.ry
                    dz = atom1.rz - atom2.rz
                    distances.append(math.sqrt(dx**2 + dy**2 + dz**2))

        return np.array(distances).reshape(dim, dim - 1)

    def inertia_tensor(self) -> NDArray[Shape["3, 3"], Float]:
        """
        Calculates the inertia tensor of the molecule.

        Returns:
            NDArray[Shape["3, 3"], Float]: inertia tensor of the molecule.
        """

        self.center_of_mass()
        Ixx = Ixy = Ixz = Iyy = Iyz = Izz = 0.0

        for atom in self.atom:
            dx = atom.rx - self.com[0]
            dy = atom.ry - self.com[1]
            dz = atom.rz - self.com[2]

            Ixx += atom.mass * (dy**2 + dz**2)
            Iyy += atom.mass * (dz**2 + dx**2)
            Izz += atom.mass * (dx**2 + dy**2)

            Ixy += atom.mass * dx * dy * -1
            Ixz += atom.mass * dx * dz * -1
            Iyz += atom.mass * dy * dz * -1

        return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

    def principal_axes(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculates the principal axes of the molecule

        Returns:
            Tuple[np.ndarray, np.ndarray]: Tuple where the first value is the Eigen Values and the second is the Eigen Vectors,
             representing the principal axes of the molecule.
        """

        try:
            evals, evecs = linalg.eigh(self.inertia_tensor())
        except ValueError:
            raise RuntimeError(
                "Error: diagonalization of inertia tensor did not converge"
            )

        return evals, evecs

    def read_position(self) -> np.ndarray:
        """Reads the position of the molecule from the position values of the atoms

        Returns:
            np.ndarray: internal position relative to atoms of the molecule
        """

        position_list = []
        for atom in self.atom:
            position_list.extend([atom.rx, atom.ry, atom.rz])
        position = np.array(position_list)
        position *= BOHR2ANG

        return position

    def update_charges(self, charges: NDArray) -> int:
        """
        Updates the charges of the atoms of the molecule and
        returns the max difference between the new and old charges
        """
        diff = 0
        for i, atom in enumerate(self.atom):
            diff = max(diff, abs(atom.chg - charges[i]))
            atom.chg = charges[i]

        return diff

    # @staticmethod
    # def update_hessian(
    #         step: np.ndarray,
    #         cur_gradient: np.ndarray,
    #         old_gradient: np.ndarray,
    #         hessian: np.ndarray,
    # ) -> np.ndarray:
    #     """
    #     Updates the Hessian of the molecule based on the current hessian, the current gradient and the previous gradient
    #
    #     Args:
    #         step (np.ndarray): step value of the iteration
    #         cur_gradient (np.ndarray): current gradient
    #         old_gradient (np.ndarray): previous gradient
    #         hessian (np.ndarray): current hessian
    #
    #     Returns:
    #         np.ndarray: updated hessian of the molecule
    #     """
    #
    #     dif_gradient = cur_gradient - old_gradient
    #
    #     mat1 = 1 / np.dot(dif_gradient, step) * np.matmul(dif_gradient.T, dif_gradient)
    #     mat2 = 1 / np.dot(step, np.matmul(hessian, step.T).T)
    #     mat2 *= np.matmul(np.matmul(hessian, step.T), np.matmul(step, hessian))
    #
    #     return hessian + mat1 - mat2

    def sizes_of_molecule(self) -> List[float]:
        """
        Calculates sides of the smallest box that the molecule could fit

        Returns:
            List[float]: list of the sizes of the molecule
        """

        x_list = []
        y_list = []
        z_list = []

        for atom in self.atom:
            x_list.append(atom.rx)
            y_list.append(atom.ry)
            z_list.append(atom.rz)

        x_max = max(x_list)
        x_min = min(x_list)
        y_max = max(y_list)
        y_min = min(y_list)
        z_max = max(z_list)
        z_min = min(z_list)

        sizes = [x_max - x_min, y_max - y_min, z_max - z_min]

        return sizes

    def standard_orientation(self) -> None:
        """
        Rotates the molecule to the standard orientation
        """

        self.center_of_mass_to_origin()
        evals, evecs = self.principal_axes()

        if round(linalg.det(evecs)) == -1:
            evecs[0, 2] *= -1
            evecs[1, 2] *= -1
            evecs[2, 2] *= -1

        if round(linalg.det(evecs)) != 1:
            raise RuntimeError(
                "Error: could not make a rotation matrix while adopting the standard orientation"
            )

        rot_matrix = evecs.T

        for atom in self.atom:
            position = np.array([atom.rx, atom.ry, atom.rz])
            new_position = np.matmul(rot_matrix, position.T).T

            atom.rx = new_position[0]
            atom.ry = new_position[1]
            atom.rz = new_position[2]

    def translate(self, vector: np.ndarray) -> "Molecule":
        """
        Creates a new Molecule object where its' atoms has been translated by a vector

        Args:
            vector (np.ndarray): translation vector

        Returns:
            Molecule: new Molecule object translated by a vector
        """

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
        inertia = self.inertia_tensor()
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

    def minimum_distance(self, molec: "Molecule") -> float:
        """
        Return the minimum distance between two molecules

        Args:
            molec (Molecule): Molecule object to be compared

        Returns:
            float: minimum distance between the two molecules
        """

        distances = []
        for atom1 in self.atom:
            if atom1.na != ghost_number:
                for atom2 in molec.atom:
                    if atom2.na != ghost_number:
                        dx = atom1.rx - atom2.rx
                        dy = atom1.ry - atom2.ry
                        dz = atom1.rz - atom2.rz
                        distances.append(math.sqrt(dx**2 + dy**2 + dz**2))

        return min(distances)
