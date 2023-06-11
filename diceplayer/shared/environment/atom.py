from diceplayer.shared.utils.ptable import atommass


class Atom:
    """
    Atom class declaration. This class is used throughout the DicePlayer program to represent atoms.

    Atributes:
        lbl (int): Dice derived variable used to represent atoms with identical energies and simetric positions.
        na (int): Atomic number of the represented atom.
        rx (float): x cartesian coordinates of the represented atom.
        ry (float): y cartesian coordinates of the represented atom.
        rz (float): z cartesian coordinates of the represented atom.
        chg (float): charge of the represented atom.
        eps (float): quantum number epsilon of the represented atom.
        sig (float): quantum number sigma of the represented atom.
    """

    def __init__(
        self,
        lbl: int,
        na: int,
        rx: float,
        ry: float,
        rz: float,
        chg: float,
        eps: float,
        sig: float,
    ) -> None:
        """
        The constructor function __init__ is used to create new instances of the Atom class.

        Args:
            lbl (int): Dice derived variable used to represent atoms with identical energies and simetric positions.
            na (int): Atomic number of the represented atom.
            rx (float): x cartesian coordinates of the represented atom.
            ry (float): y cartesian coordinates of the represented atom.
            rz (float): z cartesian coordinates of the represented atom.
            chg (float): charge of the represented atom.
            eps (float): quantum number epsilon of the represented atom.
            sig (float): quantum number sigma of the represented atom.
        """

        self.lbl = lbl
        self.na = na
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.chg = chg
        self.eps = eps
        self.sig = sig
        self.mass = atommass[self.na]
