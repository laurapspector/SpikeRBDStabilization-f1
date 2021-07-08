import io
import random
import textwrap
import unittest

import numpy as np
from dms.merge import (
    MIN_QUAL,
    bN,
    byte_array_to_str,
    str_to_byte_array,
    compare_seq_ids,
    merge_reads,
    merge_all_reads,
    read_line,
    read_seqs,
    reverse_complement,
)

def random_mismatch(base):
    return random.choice([b for b in 'ATGC' if b != base])

def random_merge_reads_test_case(min_len=150, max_len=250,
                                 min_overlap=20, max_overlap=40,
                                 min_read1_len=100, max_read1_len=120,
                                 qual_mean=25, qual_sd=5,
                                 min_mismatches=0, max_mismatches=10):
    s = ''.join(random.choices('ATGC', k=random.randint(min_len, max_len)))
    overlap = random.randint(min_overlap, max_overlap)
    s1 = s[:random.randint(min_read1_len, max_read1_len)]
    s2 = s[len(s1)-overlap:]
    q1 = np.random.normal(qual_mean, qual_sd, len(s1))\
                  .astype(np.int8)\
                  .clip(0, 40) \
                  + MIN_QUAL
    q2 = np.random.normal(qual_mean, qual_sd, len(s2))\
                  .astype(np.int8)\
                  .clip(0, 40) \
                  + MIN_QUAL
    q = np.zeros(len(s), dtype=np.int8)
    q[:len(s1)-overlap] = q1[:len(s1)-overlap]
    q[len(s1):] = q2[overlap:]
    mm_positions = sorted(random.sample(range(len(s1)-overlap, len(s1)),
                                        k=random.randint(min_mismatches,
                                                         max_mismatches)))
    s, s1, s2 = map(list, [s, s1, s2])
    for i in range(len(s1)-overlap, len(s1)):
        j = i - (len(s1) - overlap)
        if i in mm_positions:
            mm = random_mismatch(s[i])
            q[i] = abs(q1[i] - q2[j]) + MIN_QUAL
            if q1[i] > q2[j]:
                s2[j] = mm
            elif q2[j] > q1[i]:
                s1[i] = mm
            else:
                s1[i] = mm
                s[i] = 'N'
        else:
            q[i] = max(q1[i], q2[j])
    s, s1, s2 = map(lambda x: str_to_byte_array(''.join(x)), [s, s1, s2])
    s2 = reverse_complement(s2)
    q2 = np.flip(q2).copy()
    return s, q, s1, s2, q1, q2, mm_positions

def fastq_string(seqs, quals, id_suffix):
    return ''.join(
        [textwrap.dedent(
            f"""\
            @{i} {id_suffix}
            {byte_array_to_str(seq)}
            +
            {byte_array_to_str(qual)}
            """)
         for (i, (seq, qual)) in enumerate(zip(seqs, quals))])


class TestMerge(unittest.TestCase):
    def test_str_to_byte_array(self):
        self.assertTrue((str_to_byte_array('') == \
                         np.array([], dtype=np.int8)).all())
        self.assertTrue((str_to_byte_array('ATGC') == \
                         np.array([65, 84, 71, 67],
                                  dtype=np.int8)).all())
        self.assertTrue((str_to_byte_array('AATGCG') == \
                         np.array([65, 65, 84, 71, 67, 71],
                                  dtype=np.int8)).all())

    def test_byte_array_to_str(self):
        self.assertEqual(byte_array_to_str(np.array([], dtype=np.int8)), '')
        self.assertEqual(byte_array_to_str(np.array([72, 101, 108, 108, 111],
                                                    dtype=np.int8)),
                         'Hello')

    def test_reverse_complement(self):
        for s, r in [('ATGC', 'GCAT'),
                     ('TAGGACAGTA', 'TACTGTCCTA'),
                     ('', ''),
                     ('TTTTTTG', 'CAAAAAA')]:
            self.assertEqual(
                byte_array_to_str(reverse_complement(str_to_byte_array(s))),
                r)

    def test_read_line(self):
        with self.assertRaises(EOFError):
            read_line(io.StringIO(''))

        f = io.StringIO(textwrap.dedent(
            """\
            ABC 
            DEFG
            123   
            45"""))
        self.assertEqual(read_line(f), 'ABC')
        self.assertEqual(read_line(f), 'DEFG')
        self.assertEqual(read_line(f), '123')
        self.assertEqual(read_line(f), '45')
        with self.assertRaises(EOFError):
            read_line(f)

        f = io.StringIO(textwrap.dedent(
            """\
            Lorem ipsum
             dolor sit amet, 
            consectetur adipiscing 
            elit, sed do eiusmod...
            
            
            """))
        self.assertEqual(read_line(f), 'Lorem ipsum')
        self.assertEqual(read_line(f), ' dolor sit amet,')
        self.assertEqual(read_line(f), 'consectetur adipiscing')
        self.assertEqual(read_line(f), 'elit, sed do eiusmod...')
        self.assertEqual(read_line(f), '')
        self.assertEqual(read_line(f), '')
        with self.assertRaises(EOFError):
            read_line(f)

    def test_read_seqs(self):
        f = io.StringIO(textwrap.dedent(
            """\
            @1
            ATGC
            +
            AAAA
            @2
            GAATTCTTT
            +
            ABCDABCDE
            """))
        answers = [('@1', 'ATGC', '+', 'AAAA'),
                   ('@2', 'GAATTCTTT', '+', 'ABCDABCDE')]
        for res, ans in zip(read_seqs(f), answers):
            self.assertEqual(res[0], ans[0])
            self.assertTrue((res[1] == str_to_byte_array(ans[1])).all())
            self.assertEqual(res[2], ans[2])
            self.assertTrue((res[3] == str_to_byte_array(ans[3])).all())

        f = io.StringIO(textwrap.dedent(
            """\
            @longer name
            ATGCNNA
            +a
            001AABC
            @quick brown fox
            NCATTACG
            +bcd e 
            C0FFEE00
            """))
        answers = [('@longer name', 'ATGCNNA', '+a', '001AABC'),
                   ('@quick brown fox', 'NCATTACG', '+bcd e', 'C0FFEE00')]
        for res, ans in zip(read_seqs(f), answers):
            self.assertEqual(res[0], ans[0])
            self.assertTrue((res[1] == str_to_byte_array(ans[1])).all())
            self.assertEqual(res[2], ans[2])
            self.assertTrue((res[3] == str_to_byte_array(ans[3])).all())


        # Bad sequence ID.
        f = io.StringIO(textwrap.dedent(
            """\
            @good
            ATGC
            +
            BBBB
            no at
            ATGCATGC
            +
            AABBCCDD
            """))
        g = read_seqs(f)
        next(g)
        with self.assertRaises(ValueError):
            next(g)

        f = io.StringIO(textwrap.dedent(
            """\
            @s
            ATGC
            +
            BBBB
            @s2
            ATGCATGC
            +
            AABBCCD"""))
        g = read_seqs(f)
        next(g)
        with self.assertRaises(ValueError):
            next(g)

        # Bad quality ID.
        f = io.StringIO(textwrap.dedent(
            """\
            @s1
            ATGC
            +
            BBBB
            @s2
            A
            +
            A
            @s3
            ATGCATGC
            bad qual
            AABBCCDE
            """))
        g = read_seqs(f)
        next(g)
        next(g)
        with self.assertRaises(ValueError):
            next(g)

        # Bad quality ID.
        f = io.StringIO(textwrap.dedent(
            """\
            @seq 1
            ATGC
            +
            BBBB
            @seq 2
            A
             +
            A
            @seq 3
            ATGCATGC
            +
            AABBCCDE
            """))
        g = read_seqs(f)
        next(g)
        with self.assertRaises(ValueError):
            next(g)


        # Incomplete record.
        f = io.StringIO(textwrap.dedent(
            """\
            @s1
            ATGC
            +
            BBBB
            @s2
            A
            +
            A
            @s3
            ATGCATGC
            +
            """))
        g = read_seqs(f)
        next(g)
        next(g)
        with self.assertRaises(EOFError):
            next(g)

    def test_merge_reads(self):
        # Perfect matches with max qualities.
        s = str_to_byte_array('CGCGGACCTAGTCTGTAGCCGGAAGTCAAACCCAGAGTGGAGACAACATGGATTGAAAGCTTTTGACGTGCGGGGTTCGA')
        s1 = str_to_byte_array('CGCGGACCTAGTCTGTAGCCGGAAGTCAAACCCAGAGTGGAGACAACATG')
        s2 = str_to_byte_array('TCGAACCCCGCACGTCAAAAGCTTTCAATCCATGTTGTCTCCACTCTGGG')
        q = np.full(len(s1), MIN_QUAL + 40)
        sr, qr, n_mm = merge_reads(s1, s2, q, q, len(s))
        self.assertTrue((s == sr).all())
        self.assertTrue((qr == MIN_QUAL + 40).all())
        self.assertEqual(n_mm, 0)

        # Single mismatch, different read 1 and read 2 lengths.
        s = str_to_byte_array('CGGAGCCGTGAGCTCCAATCCACTTGTTATTGGAATATAC')
        s1 = str_to_byte_array('CGGAGCCGTGAGCTCCAATCCACTT')
        s2 = str_to_byte_array('GTATATTCCAATAACAAGAGGATTGGAGCT')
        q1 = np.full(len(s1), MIN_QUAL + 35)
        q1[21] = MIN_QUAL + 36
        q2 = np.full(len(s2), MIN_QUAL + 35)
        q2[18] = MIN_QUAL + 20
        q = np.full(len(s), MIN_QUAL + 35)
        q[21] = MIN_QUAL + 16
        sr, qr, n_mm = merge_reads(s1, s2, q1, q2, len(s))
        self.assertTrue((s == sr).all())
        self.assertTrue((q == qr).all())
        self.assertEqual(n_mm, 1)

        # Run many random test cases.
        for i in range(500):
            s, q, s1, s2, q1, q2, mismatches = random_merge_reads_test_case()
            sr, qr, nr = merge_reads(s1, s2, q1, q2, len(s))
            self.assertTrue((sr == s).all())
            self.assertTrue((qr == q).all())
            self.assertEqual(nr, len(mismatches))

    def test_compare_seq_ids(self):
        self.assertEqual(compare_seq_ids('@test 1',
                                         '@test 2'),
                         (True, '@test'))
        self.assertEqual(compare_seq_ids('@x:y:z 1:2:3',
                                         '@x:y:z 2:3:4'),
                         (True, '@x:y:z'))
        self.assertEqual(compare_seq_ids('@x:y:w 1:2:3',
                                         '@x:y:z 2:3:4'),
                         (False, None))
        self.assertEqual(compare_seq_ids('@x:y:z 1:2:3',
                                         '@x:y:asdf 2:3:4'),
                         (False, None))
        with self.assertRaises(ValueError):
            compare_seq_ids('@test  1', '@test 2')
        with self.assertRaises(ValueError):
            compare_seq_ids('@test 1', '@test  2')
        with self.assertRaises(ValueError):
            compare_seq_ids('@test  1', '@test  2')

    def test_merge_all_reads(self):
        # Simple test.
        f1 = textwrap.dedent(
            """\
            @first 1
            ATGCAGGACAGTGAAG
            +
            AAAAAAAAAAAAAAAA
            """)
        f2 = textwrap.dedent(
            """\
            @first 2
            ACTACTGGTCCTCACTTCACTGT
            +
            BBBBBBBBBBBBBBBBBBBBBBB
            """)
        r, *_ = merge_all_reads(io.StringIO(f1), io.StringIO(f2), 30)
        self.assertEqual(r, 'ATGCAGGACAGTGAAGTGAGGACCAGTAGT')

        # Test multiple reads.
        seqs = [
            'GACATGACTACAGGACACGAGCACTACATTTAGATACAGT',
            'AGACTAGCATCAGGCAGCACTACACAAAAGAGTGAAACCA',
            'GATCGACTAAAGCGACTACACACAGAGACTACTACGACAC',
        ]
        f1 = textwrap.dedent(
            """\
            @1 1
            GACATGACTACAGGACACGAGCACTACATT
            +
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            @2 1
            AGACTAGCATCAGGCAGCACTACAC
            +
            CCCCCCCCCCCDDDDDDDDDDDDDD
            @3 1
            GATCGACTAAAGCGACTACACACAGAGACTACT
            +
            AAAAAABBBBBABABABACCCCCCADAAACAAC
            """)
        f2 = textwrap.dedent(
            """\
            @1 2
            ACTGTATCTAAATGTAGTGCTCGTGTCCTG
            +
            999999999999999999999999999999
            @2 2
            TGGTTTCACTCTTTTGTGTAGTGCTGCCTG
            +
            BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            @3 2
            GTGTCGTAGTAGTCTCTGTG
            +
            88889999779988779999
            """)
        for s, r in zip(seqs, merge_all_reads(io.StringIO(f1),
                                              io.StringIO(f2),
                                              40)):
            self.assertEqual(s, r)

        # Test that low quality reads are dropped.
        seqs = [
            ('GACATGACTACAGGACACGAGCACTACATTTAGATACAGT', 9),
            ('AGACTAGCATCAGGCAGCACTACACAAAAGAGTGAAACCA', 13),
            ('GATCGACTAAAGCGACTACACACAGAGACTACTACGACAC', 15),
        ]
        f1 = textwrap.dedent(
            """\
            @1 1
            GACATGACTACAGGACACGAGCACTACATT
            +
            AAA*AAAAAAAAAAAAAAAAAAAAAAAAAA
            @2 1
            AGACTAGCATCAGGCAGCACTACAC
            +
            CCCCCCCCCCCDDDDDDDDDDDDDD
            @3 1
            GATCGACTAAAGCGACTACACACAGAGACTACT
            +
            AAAAAABBBBBABABABACCCCCCADAAACAAC
            """)
        f2 = textwrap.dedent(
            """\
            @1 2
            ACTGTATCTAAATGTAGTGCTCGTGTCCTG
            +
            999999999999999999999999999999
            @2 2
            TGGTTTCACTCTTTTGTGTAGTGCTGCCTG
            +
            BB.BBBBBBBBBBBBBBBBBBBBBBBBBBB
            @3 2
            GTGTCGTAGTAGTCTCTGTG
            +
            08889999779988779999
            """)
        for min_qual in [10, 13, 15, 20]:
            answers = [s for (s, q) in seqs if q >= min_qual]
            results = list(merge_all_reads(io.StringIO(f1), io.StringIO(f2), 40,
                                           min_qual=min_qual))
            self.assertEqual(len(results), len(answers))
            for ans, res in zip(answers, results):
                self.assertEqual(ans, res)

        # Test some randomly generated inputs.
        for i in range(10):
            seqs = []
            quals = []
            answers = []
            min_qual = random.randint(5, 15)
            max_mm = random.randint(0, 10)
            amplen = random.randint(150, 250)
            for i in range(1000):
                s, q, s1, s2, q1, q2, mm_positions = \
                    random_merge_reads_test_case(min_len=amplen, max_len=amplen,
                                                 min_mismatches=0,
                                                 max_mismatches=10)
                seqs.append((s1, s2))
                quals.append((q1, q2))
                if bN in s: continue
                if q.min()-MIN_QUAL < min_qual: continue
                if len(mm_positions) > max_mm: continue
                answers.append(byte_array_to_str(s))
            s1s, s2s = zip(*seqs)
            q1s, q2s = zip(*quals)
            f1 = io.StringIO(fastq_string(s1s, q1s, 1))
            f2 = io.StringIO(fastq_string(s2s, q2s, 1))
            results = list(merge_all_reads(f1, f2, amplen,
                                           max_mm=max_mm,
                                           min_qual=min_qual))
            self.assertEqual(len(results), len(answers))
            for res, ans in zip(results, answers):
                self.assertEqual(res, ans)








if __name__ == '__main__':
    unittest.main()
