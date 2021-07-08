import itertools as it

_CODONS = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '*': ['TAA', 'TAG', 'TGA']
}

_MIN_NT_DIST = {}
for aa1, aa2 in it.product(_CODONS, _CODONS):
    for c1, c2 in it.product(_CODONS[aa1], _CODONS[aa2]):
        d = sum(a != b for (a, b) in zip(c1, c2))
        _MIN_NT_DIST[aa1, aa2] = min(_MIN_NT_DIST.get((aa1, aa2), d), d)

def min_nt_dist(aa1, aa2):
    return _MIN_NT_DIST[aa1, aa2]

if __name__ == '__main__':
    examples = [('L', 'P', 1),
                ('T', 'N', 1),
                ('G', 'R', 1),
                ('V', 'P', 2),
                ('W', 'H', 3),
                ('F', '*', 2)]
    for aa1, aa2, d in examples:
        assert min_nt_dist(aa1, aa2) == min_nt_dist(aa2, aa1) == d

