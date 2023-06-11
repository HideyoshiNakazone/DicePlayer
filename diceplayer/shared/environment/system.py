from diceplayer import logger
from diceplayer.shared.environment.molecule import Molecule
from diceplayer.shared.utils.misc import BOHR2ANG
from diceplayer.shared.utils.ptable import atomsymb

import numpy as np
from numpy import linalg

import math
from copy import deepcopy
from typing import List, TextIO, Tuple


class System:
    """
    System class declaration. This class is used throughout the DicePlayer program to represent the system containing the molecules.

    Atributes:
        molecule (List[Molecule]): List of molecules of the system
        nmols (List[int]): List of number of molecules in the system
    """

    def __init__(self) -> None:
        """
        Initializes an empty system object that will be populated afterwards
        """
        self.nmols: List[int] = []
        self.molecule: List[Molecule] = []

    def add_type(self, m: Molecule) -> None:
        """
        Adds a new molecule type to the system

        Args:
            m (Molecule): The instance of the new type of molecule
        """
        if isinstance(m, Molecule) is False:
            raise TypeError("Error: molecule is not a Molecule instance")
        self.molecule.append(m)

    def update_molecule(self, position: np.ndarray) -> None:
        """Updates the position of the molecule in the Output file

        Args:
            position (np.ndarray): numpy position vector
        """

        position_in_ang = (position * BOHR2ANG).tolist()
        self.add_type(deepcopy(self.molecule[0]))

        for atom in self.molecule[-1].atom:
            atom.rx = position_in_ang.pop(0)
            atom.ry = position_in_ang.pop(0)
            atom.rz = position_in_ang.pop(0)

        rmsd, self.molecule[0] = self.rmsd_fit(-1, 0)
        self.molecule.pop(-1)

        logger.info("Projected new conformation of reference molecule with RMSD fit")
        logger.info(f"RMSD = {rmsd:>8.5f} Angstrom")

    def rmsd_fit(self, p_index: int, r_index: int) -> Tuple[float, Molecule]:
        projecting_mol = self.molecule[p_index]
        reference_mol = self.molecule[r_index]

        if len(projecting_mol.atom) != len(reference_mol.atom):
            raise RuntimeError(
                "Error in RMSD fit procedure: molecules have different number of atoms"
            )
        dim = len(projecting_mol.atom)

        new_projecting_mol = deepcopy(projecting_mol)
        new_reference_mol = deepcopy(reference_mol)

        new_projecting_mol.center_of_mass_to_origin()
        new_reference_mol.center_of_mass_to_origin()

        x = []
        y = []

        for atom in new_projecting_mol.atom:
            x.extend([atom.rx, atom.ry, atom.rz])

        for atom in new_reference_mol.atom:
            y.extend([atom.rx, atom.ry, atom.rz])

        x = np.array(x).reshape(dim, 3)
        y = np.array(y).reshape(dim, 3)

        r = np.matmul(y.T, x)
        rr = np.matmul(r.T, r)

        try:
            evals, evecs = linalg.eigh(rr)
        except:
            raise RuntimeError("Error: diagonalization of RR matrix did not converge")

        a1 = evecs[:, 2].T
        a2 = evecs[:, 1].T
        a3 = np.cross(a1, a2)

        A = np.array([a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], a3[0], a3[1], a3[2]])
        A = A.reshape(3, 3)

        b1 = np.matmul(r, a1.T).T  # or np.dot(r, a1)
        b1 /= linalg.norm(b1)
        b2 = np.matmul(r, a2.T).T  # or np.dot(r, a2)
        b2 /= linalg.norm(b2)
        b3 = np.cross(b1, b2)

        B = np.array([b1[0], b1[1], b1[2], b2[0], b2[1], b2[2], b3[0], b3[1], b3[2]])
        B = B.reshape(3, 3).T

        rot_matrix = np.matmul(B, A)
        x = np.matmul(rot_matrix, x.T).T

        rmsd = 0
        for i in range(dim):
            rmsd += (
                (x[i, 0] - y[i, 0]) ** 2
                + (x[i, 1] - y[i, 1]) ** 2
                + (x[i, 2] - y[i, 2]) ** 2
            )
        rmsd = math.sqrt(rmsd / dim)

        for i in range(dim):
            new_projecting_mol.atom[i].rx = x[i, 0]
            new_projecting_mol.atom[i].ry = x[i, 1]
            new_projecting_mol.atom[i].rz = x[i, 2]

        reference_mol.center_of_mass()

        projected_mol = new_projecting_mol.translate(reference_mol.com)

        return rmsd, projected_mol

    # def center_of_mass_distance(self, a: int, b: int) -> float:
    #     """
    #     Calculates the distance between the center of mass of two molecules
    #
    #     Args:
    #         a (Molecule): First Molecule Instance
    #         b (Molecule): Second Molecule Instance
    #
    #     Returns:
    #         float: module of the distance between the two center of masses
    #     """
    #
    #     com1 = self.molecule[a].center_of_mass()
    #     com2 = self.molecule[b].center_of_mass()
    #     dx = com1[0] - com2[0]
    #     dy = com1[1] - com2[1]
    #     dz = com1[2] - com2[2]
    #     distance = math.sqrt(dx**2 + dy**2 + dz**2)
    #
    #     return distance

    # def nearest_image(
    #     self,
    #     index_r: int,
    #     index_m: int,
    #     lx: float,
    #     ly: float,
    #     lz: float,
    #     criterium=None,
    # ) -> Tuple[float, Molecule]:
    #
    #     if criterium in None:
    #         criterium = "com"
    #
    #     if criterium != "com" and criterium != "min":
    #         raise RuntimeError("Error in value passed to function nearest_image")
    #
    #     min_dist = 1e20
    #
    #     for i in range(-1, 2):
    #         for j in range(-1, 2):
    #             for k in range(-1, 2):
    #
    #                 tr_vector = [i * lx, j * ly, k * lz]
    #                 self.add_molecule(self.molecule[index_m].translate(tr_vector))
    #
    #                 if criterium == "com":
    #                     dist = self.center_of_mass_distance(index_r, -1)
    #                 else:
    #                     dist = self.minimum_distance(index_r, -1)
    #
    #                 if dist < min_dist:
    #                     min_dist = dist
    #                     nearestmol = deepcopy(self.molecule[-1])
    #
    #     self.molecule.pop(-1)
    #
    #     return min_dist, nearestmol

    # def print_geom(self, cycle: int, fh: TextIO) -> None:
    #     """
    #             Print the geometry of the molecule in the Output file
    #
    #             Args:
    #                 cycle (int): Number of the cycle
    #                 fh (TextIO): Output file
    #             """
    #
    #     fh.write("Cycle # {}\n".format(cycle))
    #     fh.write("Number of site: {}\n".format(len(self.molecule[0].atom)))
    #     for atom in self.molecule[0].atom:
    #         symbol = atomsymb[atom.na]
    #         fh.write(
    #             "{:<2s}    {:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(
    #                 symbol, atom.rx, atom.ry, atom.rz
    #             )
    #         )
    #
    def print_charges_and_dipole(self, cycle: int) -> None:
        """
        Print the charges and dipole of the molecule in the Output file

        Args:
            cycle (int): Number of the cycle
            fh (TextIO): Output file
        """

        logger.info("Cycle # {}\n".format(cycle))
        logger.info("Number of site: {}\n".format(len(self.molecule[0].atom)))

        chargesAndDipole = self.molecule[0].charges_and_dipole()

        logger.info(
            "{:>10.6f}  {:>10.6f}  {:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(
                chargesAndDipole[0],
                chargesAndDipole[1],
                chargesAndDipole[2],
                chargesAndDipole[3],
                chargesAndDipole[4],
            )
        )
