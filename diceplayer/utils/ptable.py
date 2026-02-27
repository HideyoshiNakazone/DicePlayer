from dataclasses import dataclass
from enum import Enum


DICE_GHOST_LABEL = "Xx"

####  Number of the ghost atom
GHOST_NUMBER = 0


@dataclass(frozen=True, slots=True)
class AtomInfo:
    atomic_number: int
    symbol: str
    mass: float


class PTable(Enum):
    Xx = AtomInfo(GHOST_NUMBER, DICE_GHOST_LABEL, 0.0)
    H = AtomInfo(1, "H", 1.0079)
    He = AtomInfo(2, "He", 4.0026)
    Li = AtomInfo(3, "Li", 6.9410)
    Be = AtomInfo(4, "Be", 9.0122)
    B = AtomInfo(5, "B", 10.811)
    C = AtomInfo(6, "C", 12.011)
    N = AtomInfo(7, "N", 14.007)
    O = AtomInfo(8, "O", 15.999)
    F = AtomInfo(9, "F", 18.998)
    Ne = AtomInfo(10, "Ne", 20.180)
    Na = AtomInfo(11, "Na", 22.990)
    Mg = AtomInfo(12, "Mg", 24.305)
    Al = AtomInfo(13, "Al", 26.982)
    Si = AtomInfo(14, "Si", 28.086)
    P = AtomInfo(15, "P", 30.974)
    S = AtomInfo(16, "S", 32.065)
    Cl = AtomInfo(17, "Cl", 35.453)
    Ar = AtomInfo(18, "Ar", 39.948)
    K = AtomInfo(19, "K", 39.098)
    Ca = AtomInfo(20, "Ca", 40.078)
    Sc = AtomInfo(21, "Sc", 44.956)
    Ti = AtomInfo(22, "Ti", 47.867)
    V = AtomInfo(23, "V", 50.942)
    Cr = AtomInfo(24, "Cr", 51.996)
    Mn = AtomInfo(25, "Mn", 54.938)
    Fe = AtomInfo(26, "Fe", 55.845)
    Co = AtomInfo(27, "Co", 58.933)
    Ni = AtomInfo(28, "Ni", 58.693)
    Cu = AtomInfo(29, "Cu", 63.546)
    Zn = AtomInfo(30, "Zn", 65.409)
    Ga = AtomInfo(31, "Ga", 69.723)
    Ge = AtomInfo(32, "Ge", 72.640)
    As = AtomInfo(33, "As", 74.922)
    Se = AtomInfo(34, "Se", 78.960)
    Br = AtomInfo(35, "Br", 79.904)
    Kr = AtomInfo(36, "Kr", 83.798)
    Rb = AtomInfo(37, "Rb", 85.468)
    Sr = AtomInfo(38, "Sr", 87.620)
    Y = AtomInfo(39, "Y", 88.906)
    Zr = AtomInfo(40, "Zr", 91.224)
    Nb = AtomInfo(41, "Nb", 92.906)
    Mo = AtomInfo(42, "Mo", 95.940)
    Tc = AtomInfo(43, "Tc", 98.000)
    Ru = AtomInfo(44, "Ru", 101.07)
    Rh = AtomInfo(45, "Rh", 102.91)
    Pd = AtomInfo(46, "Pd", 106.42)
    Ag = AtomInfo(47, "Ag", 107.87)
    Cd = AtomInfo(48, "Cd", 112.41)
    In = AtomInfo(49, "In", 114.82)
    Sn = AtomInfo(50, "Sn", 118.71)
    Sb = AtomInfo(51, "Sb", 121.76)
    Te = AtomInfo(52, "Te", 127.60)
    I = AtomInfo(53, "I", 126.90)
    Xe = AtomInfo(54, "Xe", 131.29)
    Cs = AtomInfo(55, "Cs", 132.91)
    Ba = AtomInfo(56, "Ba", 137.33)
    La = AtomInfo(57, "La", 138.91)
    Ce = AtomInfo(58, "Ce", 140.12)
    Pr = AtomInfo(59, "Pr", 140.91)
    Nd = AtomInfo(60, "Nd", 144.24)
    Pm = AtomInfo(61, "Pm", 145.00)
    Sm = AtomInfo(62, "Sm", 150.36)
    Eu = AtomInfo(63, "Eu", 151.96)
    Gd = AtomInfo(64, "Gd", 157.25)
    Tb = AtomInfo(65, "Tb", 158.93)
    Dy = AtomInfo(66, "Dy", 162.50)
    Ho = AtomInfo(67, "Ho", 164.93)
    Er = AtomInfo(68, "Er", 167.26)
    Tm = AtomInfo(69, "Tm", 168.93)
    Yb = AtomInfo(70, "Yb", 173.04)
    Lu = AtomInfo(71, "Lu", 174.97)
    Hf = AtomInfo(72, "Hf", 178.49)
    Ta = AtomInfo(73, "Ta", 180.95)
    W = AtomInfo(74, "W", 183.84)
    Re = AtomInfo(75, "Re", 186.21)
    Os = AtomInfo(76, "Os", 190.23)
    Ir = AtomInfo(77, "Ir", 192.22)
    Pt = AtomInfo(78, "Pt", 195.08)
    Au = AtomInfo(79, "Au", 196.97)
    Hg = AtomInfo(80, "Hg", 200.59)
    Tl = AtomInfo(81, "Tl", 204.38)
    Pb = AtomInfo(82, "Pb", 207.20)
    Bi = AtomInfo(83, "Bi", 208.98)
    Po = AtomInfo(84, "Po", 209.00)
    At = AtomInfo(85, "At", 210.00)
    Rn = AtomInfo(86, "Rn", 222.00)
    Fr = AtomInfo(87, "Fr", 223.00)
    Ra = AtomInfo(88, "Ra", 226.00)
    Ac = AtomInfo(89, "Ac", 227.00)
    Th = AtomInfo(90, "Th", 232.04)
    Pa = AtomInfo(91, "Pa", 231.04)
    U = AtomInfo(92, "U", 238.03)
    Np = AtomInfo(93, "Np", 237.00)
    Pu = AtomInfo(94, "Pu", 244.00)
    Am = AtomInfo(95, "Am", 243.00)
    Cm = AtomInfo(96, "Cm", 247.00)
    Bk = AtomInfo(97, "Bk", 247.00)
    Cf = AtomInfo(98, "Cf", 251.00)
    Es = AtomInfo(99, "Es", 252.00)
    Fm = AtomInfo(100, "Fm", 257.00)
    Md = AtomInfo(101, "Md", 258.00)
    No = AtomInfo(102, "No", 259.00)
    Lr = AtomInfo(103, "Lr", 262.00)

    @classmethod
    def get_atomic_symbol(cls, atomic_number: int) -> str:
        for element in cls:
            if element.value.atomic_number == atomic_number:
                return element.value.symbol
        raise ValueError(f"Atomic number {atomic_number} not found in PTable.")

    @classmethod
    def get_atomic_mass(cls, atomic_number: int) -> float:
        for element in cls:
            if element.value.atomic_number == atomic_number:
                return element.value.mass
        raise ValueError(f"Atomic number {atomic_number} not found in PTable.")

    @classmethod
    def get_from_atomic_number(cls, atomic_number: int) -> AtomInfo:
        for element in cls:
            if element.value.atomic_number == atomic_number:
                return element.value
        raise ValueError(f"Atomic number {atomic_number} not found in PTable.")
