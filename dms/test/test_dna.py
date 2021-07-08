import unittest
from dms.dna import (
    is_dna_base,
    is_dna_seq,
    is_codon,
    reverse_complement,
    is_aa,
    is_aa_seq,
    translate_sequence
)

class TestDNA(unittest.TestCase):
    def test_is_dna_base(self):
        # positive examples
        self.assertTrue(is_dna_base('A'))
        self.assertTrue(is_dna_base('T'))
        self.assertTrue(is_dna_base('G'))
        self.assertTrue(is_dna_base('C'))

        # lowercase bases are not bases
        self.assertFalse(is_dna_base('a'))
        self.assertFalse(is_dna_base('g'))

        # number is not a base
        self.assertFalse(is_dna_base(0))
        self.assertFalse(is_dna_base(3.6))
        self.assertFalse(is_dna_base(-5))

        # other letters are not bases
        self.assertFalse(is_dna_base('P'))
        self.assertFalse(is_dna_base('/'))

        # other string is not base
        self.assertFalse(is_dna_base('ATGCAT'))
        self.assertFalse(is_dna_base('Pdjdf)]1i2peLL'))
        self.assertFalse(is_dna_base('!$*#//'))

        # string numbers are not bases
        self.assertFalse(is_dna_base('-1'))
        self.assertFalse(is_dna_base('100'))

        # other types are not bases
        self.assertFalse(is_dna_base(['tomorrow', 1, 4, 'hello']))
        self.assertFalse(is_dna_base({'key1':5,'key2':[4,1]}))

    def test_is_dna_seq(self):
        # positive examples
        self.assertTrue(is_dna_seq('T'))
        self.assertTrue(is_dna_seq('ATGCGATCGATATCGATATC'))
        self.assertTrue(is_dna_seq('GGCTAGGCTATTTAGCAGCTTTAGC'))

        # any other letter in a sequence is not a sequence
        self.assertFalse(is_dna_seq('ATCTAGHCTATCGAT'))
        self.assertFalse(is_dna_seq('K'))
        self.assertFalse(is_dna_seq('JIEKOFO'))

        # any sequence with a number is not a sequence
        self.assertFalse(is_dna_seq('149736287384'))
        self.assertFalse(is_dna_seq('ACTT3CGCT7ACG0CA'))

        # an integer is not a sequence
        self.assertFalse(is_dna_seq(183))
        self.assertFalse(is_dna_seq(-19))

    def test_is_codon(self):
        # three nucleotides is codon
        self.assertTrue(is_codon('ATG'))
        self.assertTrue(is_codon('TTT'))
        self.assertTrue(is_codon('CGG'))

        # other letters are not codons
        self.assertFalse(is_codon('GHA'))
        self.assertFalse(is_codon('LOG'))

        # less that or more than three is not a codon
        self.assertFalse(is_codon('A'))
        self.assertFalse(is_codon('GCTA'))
        self.assertFalse(is_codon('GAAACT'))

        # including number is not codon
        self.assertFalse(is_codon('A1A'))
        self.assertFalse(is_codon('938'))
        self.assertFalse(is_codon(1847))
        self.assertFalse(is_codon(-47))

    def test_reverse_complement(self):
        # assume all input is a sequence at this point
        # reverse compliment matches
        self.assertEqual(reverse_complement('A'), 'T')
        self.assertEqual(reverse_complement('CTGATA'), 'TATCAG')
        self.assertEqual(reverse_complement('CCGA'), 'TCGG')

        # reverse compliment does not match
        self.assertNotEqual(reverse_complement('G'), 'A')
        self.assertNotEqual(reverse_complement('ACGGA'), 'ACGGA')
        self.assertNotEqual(reverse_complement('CCATG'), 'GTACC')
        self.assertNotEqual(reverse_complement('GATT'), 'AATG')

    def test_is_aa(self):
        # it is an amino acid
        self.assertTrue(is_aa('G'))
        self.assertTrue(is_aa('Q'))

        # it is not an amino acid
        self.assertFalse(is_aa('O'))
        self.assertFalse(is_aa('U'))

        # a longer string
        self.assertFalse(is_aa('GATYE'))
        self.assertFalse(is_aa('CVAPO'))

        # numbers are not an AA
        self.assertFalse(is_aa(1947))
        self.assertFalse(is_aa('9183'))

    def test_is_aa_seq(self):
        # aa seq is aa seq
        self.assertTrue(is_aa_seq('APMNYW'))
        self.assertTrue(is_aa_seq('PL'))
        self.assertTrue(is_aa_seq('KIAQL'))

        # no aa seq is not aa seq
        self.assertFalse(is_aa_seq('OALL'))
        self.assertFalse(is_aa_seq('BXZO:'))

        # numbers return false
        self.assertFalse(is_aa_seq('194hA'))
        self.assertFalse(is_aa_seq('18331'))
        self.assertFalse(is_aa_seq(324))

    def test_translate_sequence(self):
        # a sequence is a sequence
        self.assertEqual(translate_sequence('ACGTGAACGTATACGAGAAAA'),'T*TYTRK')
        self.assertEqual(translate_sequence('GATCGCGAT'), 'DRD')

        # must be multiple of 3
        with self.assertRaises(ValueError):
            translate_sequence('AGCTA')
        with self.assertRaises(ValueError):
            translate_sequence('GATTCGGAGA')



if __name__ == '__main__':
    unittest.main()
