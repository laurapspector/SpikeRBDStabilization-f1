import unittest

from dms.tile import Tile
from dms.mutation import Mutation
from dms.main import mutation_counts

class TestMutCounts(unittest.TestCase):
    def testing_uncommon_iterables(self):
        tile = Tile(wt_seq='AGCTAC', cds_start=0, cds_end=6, first_aa=1, positions=None)
        # testing iterables that are empty, have no nts, integers...
        # empty list
        case_one = []
        # only integers
        case_two = [1,2,3,4,5,6,7,8,9,10]
        # no nucleotides
        case_three = 3*['JIPEQO'] + 4*['OPQEEW']
        # mix of nucleotides and random letters
        case_four = 4*['CTGANC'] + 8*['TTCNNA'] + 9*['JCHAYT']
        # having numbers as string
        case_five = 6*['389173'] + 9*['ACGT7G'] + ['CA56TG']
        # a weird iterable of iterables
        case_six = [8*['AGCTAC'],['CGATTA','AGGCTA','CCGGAT']]
        # an extra short or long sequence
        case_seven = 3*['ACGTGATCGAG'] + 8*['TAC'] + 99*['CGAGGAGAAA']

        self.assertEqual(mutation_counts(case_one, tile), {})

        with self.assertRaises(TypeError):
            mutation_counts(case_two, tile)
            mutation_counts(case_three, tile)
            mutation_counts(case_four, tile)
            mutation_counts(case_five, tile)
            mutation_counts(case_six, tile)
            mutation_counts(case_seven, tile)

    def testing_total_counts(self):
        tile = Tile(wt_seq='AGCGCCTAC', cds_start=0, cds_end=9, first_aa=1, positions=None)
        # testing to check that total counts are correct
        # only WT
        case_one = 83928*['AGCGCCTAC']
        case_two = 33*['AGCGCCTAC']
        # no WT
        case_three = 382*['AGCACCTAC'] + 33*['ACCGCCTAC'] + 9*['AGCGCCGAC']
        case_four = 77*['AGTTCCTAC'] + 8*['AGCGCATAC'] + 5*['AACGCCTAC'] + 193*['AGCGCCTGG']
        # WT and mutations
        case_five = 44*['AGCGCCTAC'] + 91*['ATCGCCTAC'] + 3*['AGCGCTTAC'] + 18*['AGCACCTAC'] + 99*['AGCGCCTAA']
        case_six = 991*['AGCGCCTAC'] + 183*['AGCGCCAAC'] + 13*['GGGGCCTAC'] + 38*['AGTTCCTAC'] + 21*['CGATCGAGT']

        self.assertEqual(mutation_counts(case_one, tile)[()], 83928)
        self.assertEqual(mutation_counts(case_two, tile)[()], 33)

        self.assertEqual(sum(mutation_counts(case_one, tile).values()), len(case_one))
        self.assertEqual(sum(mutation_counts(case_two, tile).values()), len(case_two))
        self.assertEqual(sum(mutation_counts(case_three, tile).values()), len(case_three))
        self.assertEqual(sum(mutation_counts(case_four, tile).values()), len(case_four))
        self.assertEqual(sum(mutation_counts(case_five, tile).values()), len(case_five))
        self.assertEqual(sum(mutation_counts(case_six, tile).values()), len(case_six))

    def testing_multiple_mutations(self):
        tile = Tile(wt_seq='GGCTATAAA', cds_start=0, cds_end=9, first_aa=1, positions=[1,2])
        # testing more than one mutation
        # two residue mutations in each and WT
        case_one = 19293*['GGCTATAAA'] + 5*['ATCTATTCC'] + 9*['GGCACGTCC'] + 58*['GCTTATCTG']
        # three residue mutations in each and no WT
        case_two = 31*['GCTCGTAAT'] + 42*['GATGACGAG'] + 77*['AATCAAGTT'] + 93*['ATGAGCACG']

        self.assertEqual(mutation_counts(case_one, tile)[()], 19293)
        with self.assertRaises(KeyError):
            mutation_counts(case_two, tile)[()]
        for i in [x for x in mutation_counts(case_one, tile).keys() if x != ()]:
            assert len(i) == 2
        for i in [x for x in mutation_counts(case_two, tile).keys() if x != ()]:
            assert len(i) == 3

    def testing_one_mutation(self):
        tile = Tile(wt_seq='AATCCCAAG', cds_start=0, cds_end=9, first_aa=1, positions=None)
        # testing some random different mutations
        # from NPK to SPK
        case_one = 8372*['AATCCCAAG'] + 918*['TCTCCCAAG']
        # from NPK to NPM
        case_two = 1173*['AATCCCAAG'] + 1113*['AATCCCATG']
        # from NPK to NQK
        case_three = 3817*['AATCCCAAG'] + 503*['AATCAAAAG']
        # from NPK to APK and to LPK
        case_four = 4839*['AATCCCAAG'] + 3892*['GCTCCCAAG'] + 1938*['CTACCCAAG']

        self.assertEqual(mutation_counts(case_one, tile)[(Mutation(pos=1, wt_aa='N', aa='S', codon='TCT'),)], 918)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='K', aa='M', codon='ATG'),)], 1113)
        self.assertEqual(mutation_counts(case_three, tile)[(Mutation(pos=2, wt_aa='P', aa='Q', codon='CAA'),)], 503)
        self.assertEqual(mutation_counts(case_four, tile)[(Mutation(pos=1, wt_aa='N', aa='A', codon='GCT'),)], 3892)
        self.assertEqual(mutation_counts(case_four, tile)[(Mutation(pos=1, wt_aa='N', aa='L', codon='CTA'),)], 1938)

    def testing_multiple_mutations(self):
        tile = Tile(wt_seq='GGGAATCGA', cds_start=0, cds_end=9, first_aa=1, positions=None)
        # testing two or three mutations
        # from GNR to ART
        case_one = 584*['GGGAATCGA'] + 113*['GCTAGAACT']
        # from GNR to YNA
        case_two = 122*['GGGAATCGA'] + 90*['TATAATGCG']
        # from GNR to GLQ
        case_three = 483*['GGGAATCGA'] + 832*['GGGTTACAG']
        # testing a mix of one mutation with multiple mutations
        # from GNR to GNA and to TYT
        case_four = 293*['GGGAATCGA'] + 910*['GGGAATGCG'] + 17*['ACTTACACC']
        # from GNR to PNR and to LNA
        case_five = 90*['GGGAATCGA'] + 18*['CCCAATCGA'] + 9*['CTCAATGCC']+ 2*['CCCAATCGA']

        self.assertEqual(len(mutation_counts(case_one, tile)), 2)
        self.assertEqual(len(mutation_counts(case_two, tile)), 2)
        self.assertEqual(len(mutation_counts(case_three, tile)), 2)
        self.assertEqual(len(mutation_counts(case_four, tile)), 3)
        self.assertEqual(len(mutation_counts(case_five, tile)), 3)
        self.assertEqual(mutation_counts(case_four, tile)[(Mutation(pos=3, wt_aa='R', aa='A', codon='GCG')),], 910)
        self.assertEqual(mutation_counts(case_five, tile)[(Mutation(pos=1, wt_aa='G', aa='P', codon='CCC')),], 20)

    def testing_different_cds_same_residue(self):
        tile = Tile(wt_seq='TATCGATCG', cds_start=0, cds_end=9, first_aa=1, positions=None)
        # testing if the same mutation with different codons does not collapse correctly
        # mutation YRS to YAS with 3 different codons
        case_one = 3712*['TATCGATCG'] + 881*['TATGCTTCG'] + 98*['TATGCGTCG'] + 3*['TATGCCTCG'] + 35*['TATAGCACG']
        # mutation YRS to YAV with 4 different codons
        case_two = 32*['TATCGATCG'] + 183*['TATCGAGTT'] + 2345*['TATCGAGTA'] + 56*['TATCGAGTC'] + 723*['TATCGAGTG']

        self.assertEqual(mutation_counts(case_one, tile)[(Mutation(pos=2, wt_aa='R', aa='A', codon='GCT')),], 881)
        self.assertEqual(mutation_counts(case_one, tile)[(Mutation(pos=2, wt_aa='R', aa='A', codon='GCG')),], 98)
        self.assertEqual(mutation_counts(case_one, tile)[(Mutation(pos=2, wt_aa='R', aa='A', codon='GCC')),], 3)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='S', aa='V', codon='GTT')),], 183)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='S', aa='V', codon='GTA')),], 2345)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='S', aa='V', codon='GTC')),], 56)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='S', aa='V', codon='GTG')),], 723)


    def testing_collapse_identical_mutations(self):
        tile = Tile(wt_seq='CGATTCGAT', cds_start=0, cds_end=9, first_aa=1, positions=None)
        # testing if many different mutations in different places in iterable collapse properly
        # total of 1 mutation, spread out, no WT
        case_one = 3*['AATTTCGAT'] + 9*['AATTTCAAA'] + ['AATTTCGAT'] + 13*['AAAAACGAT'] + 19*['AATTTCGAT']
        # total of 5 mutations, but spread out mutations across list
        case_two = 3*['CGATTGGAT'] + 4*['CGATTCAAA'] + 11*['CGATTCAAA'] + 5*['CGATTCGAT'] + 2*['CGATTCATA'] + 41*['AACTTCGAT'] + ['CGAATGGAT'] + \
                    6*['AACTTCGAT'] + 4*['CGATTCAAA'] + 9*['CGATTCATA'] + 4*['CGATTGGAT'] + ['CGAATGGAT'] + 3*['CGATTCGAT']

        self.assertEqual(mutation_counts(case_one, tile)[(Mutation(pos=1, wt_aa='R', aa='N', codon='AAT')),], 3+1+19)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='D', aa='I', codon='ATA')),], 2+9)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=3, wt_aa='D', aa='K', codon='AAA')),], 4+4+11)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=1, wt_aa='R', aa='N', codon='AAC')),], 6+41)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=2, wt_aa='F', aa='L', codon='TTG')),], 3+4)
        self.assertEqual(mutation_counts(case_two, tile)[(Mutation(pos=2, wt_aa='F', aa='M', codon='ATG')),], 1+1)
        self.assertEqual(mutation_counts(case_two, tile)[()], 5+3)


    def testing_different_start_resi(self):
        # check different start residues count correctly
        seqs = 3*['ACGGATCGA'] + 79*['AAGGATCGA'] + 103*['ACGGATTTA'] + 991*['ACGGATCGA'] + 11*['ACGGCCCGA'] + 83*['AAAGATCGA'] + \
                7*['ACGGATTTA'] + 13*['ACGGATCGA'] + 64*['GGGGATCGA'] + 21*['TTTGATCGA'] + 32*['ACGGATCCC'] + 99*['ACGAATCGA']

        tile_one = Tile(wt_seq='ACGGATCGA', cds_start=0, cds_end=9, first_aa=-1, positions=None)
        tile_two = Tile(wt_seq='ACGGATCGA', cds_start=0, cds_end=9, first_aa=0, positions=None)
        tile_three = Tile(wt_seq='ACGGATCGA', cds_start=0, cds_end=9, first_aa=53, positions=None)
        tile_four = Tile(wt_seq='ACGGATCGA', cds_start=0, cds_end=9, first_aa=910, positions=None)
        tile_five = Tile(wt_seq='ACGGATCGA', cds_start=0, cds_end=9, first_aa=423, positions=None)

        with self.assertRaises(ValueError):
            mutation_counts(seqs, tile_one)
            mutation_counts(seqs, tile_two)

        # check WT is still normal and nothing changes with new first_aa
        self.assertEqual(mutation_counts(seqs, tile_three)[()], 3+991+13)
        self.assertEqual(mutation_counts(seqs, tile_four)[()], 3+991+13)
        self.assertEqual(mutation_counts(seqs, tile_five)[()], 3+991+13)

        # check numbering is working properly
        # (chose random cases from seqs to test for each one)
        self.assertEqual(mutation_counts(seqs, tile_three)[(Mutation(pos=53, wt_aa='T', aa='K', codon='AAA')),], 83)
        self.assertEqual(mutation_counts(seqs, tile_three)[(Mutation(pos=53, wt_aa='T', aa='K', codon='AAG')),], 79)
        self.assertEqual(mutation_counts(seqs, tile_three)[(Mutation(pos=55, wt_aa='R', aa='L', codon='TTA')),], 103+7)
        self.assertEqual(mutation_counts(seqs, tile_four)[(Mutation(pos=910, wt_aa='T', aa='G', codon='GGG')),], 64)
        self.assertEqual(mutation_counts(seqs, tile_four)[(Mutation(pos=911, wt_aa='D', aa='A', codon='GCC')),], 11)
        self.assertEqual(mutation_counts(seqs, tile_five)[(Mutation(pos=424, wt_aa='D', aa='N', codon='AAT')),], 99)
        self.assertEqual(mutation_counts(seqs, tile_five)[(Mutation(pos=425, wt_aa='R', aa='P', codon='CCC')),], 32)

    def testing_frameshift(self):
        # check frameshift works properly
        seqs = 55*['ACGATCGATC'] + 91*['ACGCTCGATC'] + 12*['ACGATCGATG']

        tile_one = Tile(wt_seq='ACGATCGATC', cds_start=1, cds_end=10, first_aa=1, positions=None)
        tile_two = Tile(wt_seq='ACGATCGATC', cds_start=0, cds_end=9, first_aa=1, positions=None)
        tile_three = Tile(wt_seq='ACGATCGATC', cds_start=2, cds_end=8, first_aa=1, positions=None)
        tile_four = Tile(wt_seq='ACGATCGATC', cds_start=3, cds_end=9, first_aa=1, positions=None)
        tile_five = Tile(wt_seq='ACGATCGATC', cds_start=4, cds_end=10, first_aa=1, positions=None)

        # check WT is correct
        self.assertEqual(mutation_counts(seqs, tile_one)[()], 55)
        self.assertEqual(mutation_counts(seqs, tile_two)[()], 55)
        self.assertEqual(mutation_counts(seqs, tile_three)[()], 55)
        self.assertEqual(mutation_counts(seqs, tile_four)[()], 55)
        self.assertEqual(mutation_counts(seqs, tile_five)[()], 55)

        # check mutations
        # randomly chosen from above seqs list
        self.assertEqual(mutation_counts(seqs, tile_two)[(Mutation(pos=2, wt_aa='I', aa='L', codon='CTC')),], 91)
        self.assertEqual(mutation_counts(seqs, tile_three)[(Mutation(pos=1, wt_aa='D', aa='A', codon='GCT')),], 91)
        self.assertEqual(mutation_counts(seqs, tile_four)[(Mutation(pos=1, wt_aa='I', aa='L', codon='CTC')),], 91)
        self.assertEqual(mutation_counts(seqs, tile_five)[(Mutation(pos=2, wt_aa='I', aa='M', codon='ATG')),], 12)

if __name__ == '__main__':
    unittest.main()
