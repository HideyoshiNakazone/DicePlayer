from diceplayer.environment.molecule import Molecule

from typing_extensions import List

from dataclasses import dataclass, field


@dataclass
class System:
    """
    System class declaration. This class is used throughout the DicePlayer program to represent the system containing the molecules.

    Atributes:
        molecule (List[Molecule]): List of molecules of the system
        nmols (List[int]): List of number of molecules in the system
    """

    nmols: List[int] = field(default_factory=list)
    molecule: List[Molecule] = field(default_factory=list)

    def add_type(self, m: Molecule) -> None:
        """
        Adds a new molecule type to the system

        Args:
            m (Molecule): The instance of the new type of molecule
        """
        if not isinstance(m, Molecule):
            raise TypeError("Error: molecule is not a Molecule instance")
        self.molecule.append(m)
