import unittest
from dms.mutation import (
    Mutation,
    NontargetMutation,
    WildType,
    AminoAcidMutation
)

class TestMutation(unittest.TestCase):
    def test_construction(self):
        # pos must be an integer.
        with self.assertRaises(TypeError):
            Mutation(10.0, 'R', 'H', 'CAT')
        with self.assertRaises(TypeError):
            Mutation('thirty', 'I', 'K', 'AAA')

        # pos must be positive.
        with self.assertRaises(ValueError):
            Mutation(-100, 'G', 'Y', 'TAC')
        with self.assertRaises(ValueError):
            Mutation(0, 'T', 'S', 'AGT')

        # wt_aa must be an AA.
        with self.assertRaises(TypeError):
            Mutation(35, 'X', 'G', 'GGA')
        with self.assertRaises(TypeError):
            Mutation(42, 'ATG', 'A', 'GCT')

        # aa must be an AA.
        with self.assertRaises(TypeError):
            Mutation(5, 'Y', 12, 'GCA')
        with self.assertRaises(TypeError):
            Mutation(98, 'A', 'dog', 'GCT')

        # codon must be a codon.
        with self.assertRaises(TypeError):
            Mutation(283, 'D', 'E', 'UUG')
        with self.assertRaises(TypeError):
            Mutation(134, 'G', 'K', 79)

        # codon must match aa.
        with self.assertRaises(ValueError):
            Mutation(31, 'C', 'T', 'GGG')
        with self.assertRaises(ValueError):
            Mutation(1, 'W', 'K', 'ATG')

    def test_repr(self):
        self.assertEqual(repr(Mutation(423, 'M', 'I', 'ATC')), 'M423I-ATC')
        self.assertEqual(repr(Mutation(931, 'V', 'K', 'AAG')), 'V931K-AAG')


class TestNontargetMutation(unittest.TestCase):
    def test_construction(self):
        # pos must be an integer.
        with self.assertRaises(TypeError):
            NontargetMutation(1.0, 'A', 'C')
        with self.assertRaises(TypeError):
            NontargetMutation('five', 'T', 'A')

        # pos must be non-negative.
        with self.assertRaises(ValueError):
            NontargetMutation(-10, 'G', 'C')
        with self.assertRaises(ValueError):
            NontargetMutation(-2041, 'A', 'G')

        # wt_base must be a base.
        with self.assertRaises(TypeError):
            NontargetMutation(22, 'X', 'T')
        with self.assertRaises(TypeError):
            NontargetMutation(1032, 5, 'G')

        # base must be a base.
        with self.assertRaises(TypeError):
            NontargetMutation(43, 'G', 'CAT')
        with self.assertRaises(TypeError):
            NontargetMutation(935, 'T', None)

        # wt_base and base cannot be equal.
        with self.assertRaises(ValueError):
            NontargetMutation(28, 'A', 'A')
            NontargetMutation(256, 'C', 'C')

class TestAminoAcidMutation(unittest.TestCase):
    def test_construction(self):
        # pos must be an integer.
        with self.assertRaises(TypeError):
            AminoAcidMutation(-0.5, 'Y', 'W')
        with self.assertRaises(TypeError):
            AminoAcidMutation('twenty', 'C', 'G')

        # pos must be positive.
        with self.assertRaises(ValueError):
            AminoAcidMutation(-20, 'G', 'E')
        with self.assertRaises(ValueError):
            AminoAcidMutation(0, 'H', 'K')

        # wt_aa must be an AA.
        with self.assertRaises(TypeError):
            AminoAcidMutation(39, 9, 'E')
        with self.assertRaises(TypeError):
            AminoAcidMutation(24, 'X', 'G')

        # aa must be an AA.
        with self.assertRaises(TypeError):
            AminoAcidMutation(722, 'G', 'B')
        with self.assertRaises(TypeError):
            AminoAcidMutation(4, 'L', None)

    def test_repr(self):
        self.assertEqual(repr(AminoAcidMutation(48, 'D', 'K')), 'D48K')
        self.assertEqual(repr(AminoAcidMutation(291, 'G', 'H')), 'G291H')

    def test_from_mutation(self):
        m = Mutation(73, 'I', 'L', 'TTG')
        aam = AminoAcidMutation.from_mutation(m)
        self.assertEqual((aam.pos, aam.wt_aa, aam.aa), (73, 'I', 'L'))
        m = Mutation(228, 'M', 'K', 'AAA')
        aam = AminoAcidMutation.from_mutation(m)
        self.assertEqual((aam.pos, aam.wt_aa, aam.aa), (228, 'M', 'K'))


if __name__ == '__main__':
    unittest.main()
