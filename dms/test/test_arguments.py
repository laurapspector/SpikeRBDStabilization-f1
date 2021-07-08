import unittest
import mock
import io
import textwrap
import argparse
import configparser

from dms.tile import Tile
from dms.arguments import (ARGUMENTS,
                       parse_params,
                       parse_tiles,
                       parse_samples,
                       parse_experiments,
                       parse_proteins,
                       parse_config,
                       parse_args_and_read_config)

class TestArguments(unittest.TestCase):
    @mock.patch('configparser.open')
    def testing_parse_params(self, mockFileOpen: mock.MagicMock):
        # checking normal parameters section
        f1 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: 10
        min_quality: 10
        fastq_file_dir: 'testpath'
        """
        ))

        f2 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: 0
        min_quality: 20
        fastq_file_dir: 'new_path'
        """
        ))

        # check error for string for max mismatches
        f3 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: '10'
        min_quality: 10
        fastq_file_dir: 'testpath'
        """
        ))

        # check error for int with file directory
        f4 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: 10
        min_quality: 10
        fastq_file_dir: 10
        """
        ))

        # check error for adding a new unknown option
        f5 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: 10
        min_quality: 10
        fastq_file_dir: 'testpath'
        new_option: 8
        """
        ))

        # if no params section
        f6 = io.StringIO(textwrap.dedent(
        """\
        [Random Section]
        max_mismatches: 10
        min_quality: 10
        fastq_file_dir: 'testpath'
        """
        ))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f1)
        self.assertEqual(parse_params(ARGUMENTS, config), argparse.Namespace(fastq_file_dir='testpath', max_mismatches=10, min_quality=10))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f2)
        self.assertEqual(parse_params(ARGUMENTS, config), argparse.Namespace(fastq_file_dir='new_path', max_mismatches=0, min_quality=20))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f6)
        self.assertEqual(parse_params(ARGUMENTS, config), argparse.Namespace())

        with self.assertRaises(argparse.ArgumentTypeError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f3)
            parse_params(ARGUMENTS, config)

            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f4)
            parse_params(ARGUMENTS, config)

        with self.assertRaises(ValueError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f5)
            parse_params(ARGUMENTS, config)

    @mock.patch('configparser.open')
    def testing_parse_tiles(self, mockFileOpen: mock.MagicMock):
        # check that normal tile will work
        f1 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'TGGAGGAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 3
        """
        ))

        # a different normal tile
        f2 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T2]
        wt_seq: 'CGATCGATT'
        first_aa: 1
        cds_start: 2
        cds_end: 8
        positions: 1, 2
        """
        ))

        # check no error if WT is longer than cds_end
        f3 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T3]
        wt_seq: 'GGATCGAAT'
        first_aa: 1
        cds_start: 0
        cds_end: 6
        positions: 1, 2
        """
        ))

        # check error if cds_start and cds_end are not multiple of three
        f4 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'CGATTTTAC'
        first_aa: 2
        cds_start: 1
        cds_end: 5
        positions: 1, 2
        """
        ))

        # check error if WT is shorter than cds_end-cds_start
        f5 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T3]
        wt_seq: 'GAAGGC'
        first_aa: 1
        cds_start: 1
        cds_end: 10
        positions: 1, 2
        """
        ))

        # check error if more positions than possible
        f6 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T8]
        wt_seq: 'CTATGC'
        first_aa: 1
        cds_start: 0
        cds_end: 6
        positions: 1, 2, 3
        """
        ))

        # check error for no tile section
        f7 = io.StringIO(textwrap.dedent(
        """\

        """
        ))

        # check error for empty sections tiles
        f8 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T2]
        wt_seq: 'GGATCG'
        first_aa: 1
        cds_start:
        cds_end: 6
        positions: 1, 2
        """
        ))

        # test multiple tiles in one file
        f9 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'AGCTTC'
        first_aa: 19
        cds_start: 0
        cds_end: 6
        positions: 19, 20

        [Tile:T4]
        wt_seq: 'GGATCG'
        first_aa: 2
        cds_start: 0
        cds_end: 6
        positions: 2, 3

        [Tile:T7]
        wt_seq: 'AGCTTCGAT'
        first_aa: 8
        cds_start: 0
        cds_end: 9
        positions: 8, 10
        """
        ))

        # checking no error for adding random new option
        f10 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'CGTAGCTAG'
        first_aa: 1
        cds_start: 0
        new_option: 10
        cds_end: 9
        positions: 1, 2
        """
        ))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f1)
        self.assertEqual(parse_tiles(config)['T1'], Tile('TGGAGGAGG',1,0,9,[1,3]))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f2)
        self.assertEqual(parse_tiles(config)['T2'], Tile('CGATCGATT',1,2,8,[1,2]))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f3)
        self.assertEqual(parse_tiles(config)['T3'], Tile('GGATCGAAT',1,0,6,[1,2]))

        with self.assertRaises(ValueError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f4)
            parse_tiles(config)

            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f5)
            parse_tiles(config)

            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f6)
            parse_tiles(config)

        with self.assertRaises(argparse.ArgumentTypeError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f7)
            parse_tiles(config)

        with self.assertRaises(ValueError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f8)
            parse_tiles(config)

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f9)
        self.assertEqual(list(parse_tiles(config).keys()), ['T1', 'T4', 'T7'])

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f10)
        self.assertEqual(parse_tiles(config)['T1'], Tile('CGTAGCTAG',1,0,9,[1,2]))

    @mock.patch('configparser.open')
    def testing_parse_samples(self, mockFileOpen: mock.MagicMock):
        # checking normal samples section
        f1 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'TGGAGGAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 3

        [Tile:T2]
        wt_seq: 'AAGCTA'
        first_aa: 1
        cds_start: 0
        cds_end: 6
        positions: 1, 2

        [Samples]
        1_A1:    'T1', 'Sample_File1.fastq.gz', 'Sample_File2.fastq.gz'
        1_Display:    'T1', 'Sample_File3.fastq.gz', 'Sample_File4.fastq.gz'
        1_Control: 'T1', 'Sample_File5.fastq.gz', 'Sample_File6.fastq.gz'
        2_A1:    'T2', 'Sample_File7.fastq.gz', 'Sample_File8.fastq.gz'
        2_Display:    'T2', 'Sample_File9.fastq.gz', 'Sample_File10.fastq.gz'
        2_Control: 'T2', 'Sample_File11.fastq.gz', 'Sample_File12.fastq.gz'
        """
        ))

        # checking another normal samples section
        f2 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Tile:T2]
        wt_seq: 'GCTAGCTAT'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Samples]
        1_Display:    'T1', 'Sample_File3.fastq.gz', 'Sample_File4.fastq.gz'
        1_Control: 'T1', 'Sample_File5.fastq.gz', 'Sample_File6.fastq.gz'
        2_Display:    'T2', 'Sample_File9.fastq.gz', 'Sample_File10.fastq.gz'
        2_Control: 'T2', 'Sample_File11.fastq.gz', 'Sample_File12.fastq.gz'
        """
        ))

        # check error if tile specified is different that T in sample
        f3 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Tile:T4]
        wt_seq: 'GCTAGCTAT'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Samples]
        1_B2:    'T1', 'Sample_File3.fastq.gz', 'Sample_File4.fastq.gz'
        1_C2: 'T1', 'Sample_File5.fastq.gz', 'Sample_File6.fastq.gz'
        2_C2:    'T2', 'Sample_File9.fastq.gz', 'Sample_File10.fastq.gz'
        2_D2: 'T2', 'Sample_File11.fastq.gz', 'Sample_File12.fastq.gz'
        """
        ))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f1)
        tiles = parse_tiles(config)
        self.assertEqual(parse_samples(config, tiles)['2_A1'], ('T2', ('Sample_File7.fastq.gz', 'Sample_File8.fastq.gz')))
        self.assertEqual(parse_samples(config, tiles)['1_Display'], ('T1', ('Sample_File3.fastq.gz', 'Sample_File4.fastq.gz')))
        self.assertEqual(parse_samples(config, tiles)['2_Control'], ('T2', ('Sample_File11.fastq.gz', 'Sample_File12.fastq.gz')))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f2)
        tiles = parse_tiles(config)
        self.assertEqual(parse_samples(config, tiles)['1_Control'], ('T1', ('Sample_File5.fastq.gz', 'Sample_File6.fastq.gz')))
        self.assertEqual(parse_samples(config, tiles)['1_Display'], ('T1', ('Sample_File3.fastq.gz', 'Sample_File4.fastq.gz')))
        self.assertEqual(parse_samples(config, tiles)['2_Control'], ('T2', ('Sample_File11.fastq.gz', 'Sample_File12.fastq.gz')))

        with self.assertRaises(ValueError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f3)
            tiles = parse_tiles(config)
            print(parse_samples(config, tiles))


    @mock.patch('configparser.open')
    def testing_parse_experiments(self, mockFileOpen: mock.MagicMock):
        # no experiments section
        f1 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: 10
        min_quality: 10
        fastq_file_dir: 'testpath'
        """
        ))

        # make sure normal config works well
        f2 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        1_12-1: '1_Ref', '1_12-1'
        """
        ))

        # any experiment name is acceptable
        f3 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        Any_Name1: '1_Ref', '1_6-29'
        Any_Name2: '1_Ref', '1_12-1'
        """
        ))

        # if experiment not in samples raise error
        f4 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        1_12-1: '1_Ref', '1_12-4'
        """
        ))

        # if experiments uses files from different tiles
        f4 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'
        2_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        2_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        2_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '2_Ref', '1_6-29'
        1_12-1: '1_Ref', '1_12-4'
        """
        ))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f1)
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_experiments(config, ('T2', ('Sample_File7.fastq.gz', 'Sample_File8.fastq.gz')))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f2)
        self.assertEqual(parse_experiments(config, parse_samples(config,parse_tiles(config))), {'1_6-29': ('1_Ref', '1_6-29'), '1_12-1': ('1_Ref', '1_12-1')})

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f3)
        self.assertEqual(parse_experiments(config, parse_samples(config,parse_tiles(config))), {'Any_Name1': ('1_Ref', '1_6-29'), 'Any_Name2': ('1_Ref', '1_12-1')})

        with self.assertRaises(ValueError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f4)
            parse_experiments(config, parse_samples(config,parse_tiles(config)))

            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f5)
            parse_experiments(config, parse_samples(config,parse_tiles(config)))


    @mock.patch('configparser.open')
    def testing_parse_proteins(self, mockFileOpen: mock.MagicMock):
        # no proteins section
        f1 = io.StringIO(textwrap.dedent(
        """\
        [Parameters]
        max_mismatches: 10
        min_quality: 10
        fastq_file_dir: 'testpath'
        """
        ))

        # empty proteins raises error
        f2 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2, 3

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        1_12-1: '1_Ref', '1_12-1'

        [Proteins]

        """
        ))

        # make sure normal one works
        f3 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Tile:T2]
        wt_seq: 'GCTAGCTAGAAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 12
        positions: 3, 4

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'
        2_6-29:    'T2', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        2_12-1:    'T2', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        2_Ref:     'T2', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        2_6-29: '2_Ref', '2_6-29'
        1_12-1: '1_Ref', '1_12-1'

        [Proteins]
        6-29: '1_6-29', '2_6-29'
        """
        ))

        # experiment is not defined error check
        f4 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Tile:T2]
        wt_seq: 'GCTAGCTAGAAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 12
        positions: 3, 4

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'
        2_6-29:    'T2', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        2_12-1:    'T2', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        2_Ref:     'T2', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        1_12-1: '1_Ref', '1_12-1'

        [Proteins]
        6-29: '1_6-29', '2_6-29'
        """
        ))

        # cross over for positions should raise an error
        f5 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Tile:T2]
        wt_seq: 'GCTAGCTAGAAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 12
        positions: 2, 3

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'
        2_6-29:    'T2', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        2_12-1:    'T2', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        2_Ref:     'T2', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        2_6-29: '2_Ref', '2_6-29'
        1_12-1: '1_Ref', '1_12-1'

        [Proteins]
        6-29: '1_6-29', '2_6-29'
        """
        ))

        # random name does not matter for experiment or proteins
        f7 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Tile:T2]
        wt_seq: 'GCTAGCTAGAAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 12
        positions: 3, 4

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'
        2_6-29:    'T2', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        2_12-1:    'T2', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        2_Ref:     'T2', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        Random1: '1_Ref', '1_6-29'
        Random2: '2_Ref', '2_6-29'
        1_12-1: '1_Ref', '1_12-1'

        [Proteins]
        mAb: 'Random1', 'Random2'
        """
        ))

        # multiple proteins work
        f8 = io.StringIO(textwrap.dedent(
        """\
        [Tile:T1]
        wt_seq: 'GCTAGCTAGA'
        first_aa: 1
        cds_start: 0
        cds_end: 9
        positions: 1, 2

        [Tile:T2]
        wt_seq: 'GCTAGCTAGAAGG'
        first_aa: 1
        cds_start: 0
        cds_end: 12
        positions: 3, 4

        [Samples]
        1_6-29:    'T1', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        1_12-1:    'T1', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        1_Ref:     'T1', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'
        2_6-29:    'T2', '1-6-29_S6_L001_R1_001.fastq.gz', '1-6-29_S6_L001_R2_001.fastq.gz'
        2_12-1:    'T2', '1-12-1_S4_L001_R1_001.fastq.gz', '1-12-1_S4_L001_R2_001.fastq.gz'
        2_Ref:     'T2', '1U_S10_L001_R1_001.fastq.gz',    '1U_S10_L001_R2_001.fastq.gz'

        [Experiments]
        1_6-29: '1_Ref', '1_6-29'
        2_6-29: '2_Ref', '2_6-29'
        1_12-1: '1_Ref', '1_12-1'
        2_12-1: '2_Ref', '2_12-1'

        [Proteins]
        6-29: '1_6-29', '2_6-29'
        12-1: '1_12-1', '2_12-1'
        """
        ))

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f1)
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_proteins(config, ('T2', ('Sample_File7.fastq.gz', 'Sample_File8.fastq.gz')), Tile('CGTAGCTAG',1,0,9,[1,2]), {'1_6-29': ('1_Ref', '1_6-29'), '1_12-1': ('1_Ref', '1_12-1')})


        with self.assertRaises(ValueError):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f2)
            parse_proteins(parse_proteins(config, parse_tiles(config), parse_samples(config, parse_tiles(config)), parse_experiments(config, parse_samples(config, parse_tiles(config)))))

            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f4)
            parse_proteins(parse_proteins(config, parse_tiles(config), parse_samples(config, parse_tiles(config)), parse_experiments(config, parse_samples(config, parse_tiles(config)))))

            config = configparser.ConfigParser()
            config.optionxform = str
            config.read_file(f5)
            parse_proteins(parse_proteins(config, parse_tiles(config), parse_samples(config, parse_tiles(config)), parse_experiments(config, parse_samples(config, parse_tiles(config)))))


        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f3)
        self.assertEqual(parse_proteins(config, parse_tiles(config), parse_samples(config, parse_tiles(config)), parse_experiments(config, parse_samples(config, parse_tiles(config)))), {'6-29': ('1_6-29', '2_6-29')})

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f7)
        self.assertEqual(parse_proteins(config, parse_tiles(config), parse_samples(config, parse_tiles(config)), parse_experiments(config, parse_samples(config, parse_tiles(config)))), {'mAb': ('Random1', 'Random2')})

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read_file(f8)
        self.assertEqual(parse_proteins(config, parse_tiles(config), parse_samples(config, parse_tiles(config)), parse_experiments(config, parse_samples(config, parse_tiles(config)))), {'6-29': ('1_6-29', '2_6-29'), '12-1': ('1_12-1', '2_12-1')})


if __name__ == '__main__':
    unittest.main()
