from diceplayer.shared.environment.atom import Atom

import unittest


class TestAtom(unittest.TestCase):
    def test_class_instantiation(self):
        atom = Atom(
            lbl=1,
            na=1,
            rx=1.0,
            ry=1.0,
            rz=1.0,
            chg=1.0,
            eps=1.0,
            sig=1.0,
        )

        self.assertIsInstance(atom, Atom)