import sys
from six_frame_translation import six_frame_translation
from collections import defaultdict
from utils import *
import numpy as np

KMER_FILE = 'kmer_dict_k_4_num_prots_2.pickle'
UNIPROT_FILE = 'uniprot_sprot.fasta'

def trace(dp, x, y)->int:
    '''
    Backtrace edit-distance matrix D for strings x and y
    Returns the index at which the trace started

    Inputs
    dp -> list[list]
    x, y -> strings

    Output:
    ind -> int
    '''

    i, j = len(x), len(y)
    while i > 0:
        diag, vert, horz = float('inf'), float('inf'), float('inf')
        delta = None
        if i > 0 and j > 0:
            delta = 0 if x[i-1] == y[j-1] else 1
            diag = dp[i-1, j-1] + delta
        if i > 0:
            vert = dp[i-1, j] + 1
        if j > 0:
            horz = dp[i, j-1] + 1
        if diag <= vert and diag <= horz:
            # diagonal was best
            i -= 1; j -= 1
        elif vert <= horz:
            # vertical was best; this is an insertion in x w/r/t y
            i -= 1
        else:
            # horizontal was best
            j -= 1

    # j is index of the first (leftmost) character of t involved in the
    # alignment
    return j


def dp_edit(p,t):
    '''
    Returns edit distance between pattern p and string t
    Will return left most match if there are multiple
    '''
    dp = np.zeros((len(p)+1, len(t)+1), dtype = int)
    dp[1:,0] = range(1, len(p)+1)


    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            delta = 1 if p[i-1] != t[j-1] else 0
            dp[i, j] = min(dp[i-1, j-1] + delta, dp[i-1, j] + 1, dp[i, j-1] + 1)

    # get minimum edit distance and it's col number in last row

    min_ed_index = np.argmin(dp[-1,:])
    min_ed = dp[-1, min_ed_index]

    trace_start_ind = trace(dp, p, t[:min_ed_index])
    return trace_start_ind, min_ed


def query(read, max_mismatch = 4):
    """ Function that outputs (approximate) matches of a DNA read against
    a protein database. Implemented using k-mers and the pigeonhole principle.

    Assumes max number of mismatches is fixed (4)
    Assumes length of read is fixed (60)

    Does not account for gaps or insertions (only mismatches)
    """
    assert (len(read) == 60)

    # kmer_dict holds mapping from {kmer: [(protein_id, offset)]}
    kmer_dict = read_from_pickle(KMER_FILE)

    # protein_dict holds mapping from {protein_id : protein_seq}
    protein_dict = create_protein_map(UNIPROT_FILE, 1)

    # get amino acid sequence by six frame translation
    prot_translations = six_frame_translation(read)
    print(prot_translations)
    # Result stored as a dict of tuples (id, offset)->edit_distance
    occurrences = {}

    for protein in prot_translations:
        for i in range(len(protein)//4): # for all 15 partitions

            # for each partition
            part_off = 4*i
            part = protein[part_off:part_off+4]

            if part not in kmer_dict:
                continue

            # hits is a list of [(protein_id, offset_in_protein)]
            hits = kmer_dict[part]


            for (id, offset) in hits:
                t = protein_dict[id]

                if len(t) == 0:
                    continue

                '''count_mismatch = 0

                if offset - part_off + len(protein) <= len(t) and offset - part_off >= 0:
                    for j in range(0, len(protein)):
                        if t[offset-part_off+j] != protein[j]:
                            count_mismatch += 1
                            if count_mismatch > max_mismatch:
                                break
                    if count_mismatch <= max_mismatch:
                        occurrences.add((id, offset-part_off))
                        break # now check the next possible protein in database
                '''

                # left edge of T to include in DP matrix
                lf = max(0, offset - part_off - max_mismatch)

                # right edge of T to include in DP matrix
                rt = min(len(t), offset - part_off + len(protein) + max_mismatch)

                # get min edit distance, and start offset of traceback
                start_off, min_ed = dp_edit(protein, t[lf:rt])

                start_off += lf # get actual offset in original protein

                if min_ed <= max_mismatch:
                    occurrences[(id, start_off)] = min_ed


    return occurrences

# Test Cases
def main():

    # test protein = "ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG"
    # pid = "Q6GZX4"

    # 1 Mismatch - should pass
    occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG")
    print(occ)
    assert(("Q6GZX4", 0) in occ)

    # 2 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG")
    assert(("Q6GZX4", 0) in occ)

    # 3 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTTCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG")
    assert(("Q6GZX4", 0) in occ)

    # 4 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTTCTGAGAGAATATGATCGCCGCCGCCGCATGGAAGCG")
    assert(("Q6GZX4", 0) in occ)

    # > 4 Mismatches - should fail
    occ = query("TTGCTTTAGCGCGGAAGATGTTCTGAGAGAAGTATGATCGCCGCCGCCACATGGAAGCTG")
    assert(("Q6GZX4", 0) not in occ)


if __name__ == "__main__":
    main()
