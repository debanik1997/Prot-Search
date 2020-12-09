'''
This python code matches given DNA read with protein database using a kmer index
and a dp based edit distance

Written for Group 63 Project, Computational Genomics Fa 2020

Usage:

python index_assisted_dp.py
'''

import sys
from six_frame_translation import six_frame_translation
from collections import defaultdict
import numpy as np
import pickle

KMER_DICT_FILE = '../../data/kmer_dict_k_4_num_prots_100.pickle'
PROTEIN_DICT_FILE = '../../data/protein_dict_num_prots_100.pickle'

def read_from_pickle(f):
    """ Reads a dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        dictionary = pickle.load(fd)
    return dictionary

def trace(dp, x, t)->int:
    '''
    Backtrace edit-distance matrix D for strings x and y
    Returns the index at which the trace started

    Inputs
    dp -> list[list]
    x, y -> strings

    Output:
    ind -> int
    '''

    i, j = len(x), len(t)
    while i > 0:
        diag, vert, horz = float('inf'), float('inf'), float('inf')
        delta = None
        if i > 0 and j > 0:
            delta = 0 if x[i-1] == t[j-1] else 1
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


def query(read, kmer_dict, protein_dict, max_mismatch = 4, l = 4):
    """ Function that outputs (approximate) matches of a DNA read against
    a protein database. Implemented using k-mers and the pigeonhole principle.

    Default value of l-mer length is 4
    Default value of max mismatches allowed is 4

    """

    #will throw AssertionError if we can't make >max_mismatch partitions of size l
    '''
    Min # of partitions needed for pigeonhole = (l*max_mismatch+1)
    Min read length = min # partitions * 3
    '''
    min_rl = 3*(l*max_mismatch + 1)
    assert len(read)>= min_rl, f'Please enter a read of length {min_rl} or higher'

    # get amino acid sequence by six frame translation
    prot_translations = six_frame_translation(read)
    # Result stored as a dict of tuples (id, offset)->edit_distance
    occurrences = {}
    seen = set()

    for protein in prot_translations:
        for i in range(len(protein) - l + 1): # for all 5 partitions

            # for each partition
            part_off = i
            part = protein[part_off:part_off+l]

            if part not in kmer_dict:
                continue

            # hits is a list of [(protein_id, offset_in_protein)]
            hits = kmer_dict[part]

            for (id, offset) in hits:
                t = protein_dict[id]

                if len(t) == 0:
                    continue

                '''
                Old code: just here for reference
                count_mismatch = 0

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

                # check if we have seen this substring before
                if (id, lf, rt) in seen:
                    continue

                seen.add((id, lf, rt))

                # get min edit distance, and start offset of traceback
                start_off, min_ed = dp_edit(protein, t[lf:rt])
                start_off += lf # get actual offset in original protein

                if min_ed <= max_mismatch:
                    occurrences[(id, start_off)] = min_ed

    return occurrences

# Test Cases
def main():

    # kmer_dict holds mapping from {kmer: [(protein_id, offset)]}
    kmer_dict = read_from_pickle(KMER_DICT_FILE)

    # protein_dict holds mapping from {protein_id : protein_seq}
    protein_dict = read_from_pickle(PROTEIN_DICT_FILE)


    #########################################
    # Testing on some inputs
    #########################################
    # test protein = "ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG"
    # pid = "Q6GZX4"

    # length = 50 - should fail
    try:
        occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCG", kmer_dict, protein_dict)
    except AssertionError as e:
        # will print assertion error
        print(e, '. Continuing with other tests.')

    # 1 insertion and 1 deletion - should pass
    occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGCG", kmer_dict, protein_dict)
    assert ("Q6GZX4", 0) in occ, print(occ)

    # 1 Mismatch - should pass
    occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict)
    assert(("Q6GZX4", 0) in occ)

    # 2 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict)
    assert(("Q6GZX4", 0) in occ)

    # 3 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTTCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict)
    assert(("Q6GZX4", 0) in occ)

    # 4 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTTCTGAGAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict)
    assert(("Q6GZX4", 0) in occ)

    # > 4 Mismatches - should fail
    occ = query("TTGCTTTAGCGCGGAAGATGTTCTGAGAGAAGTATGATCGCCGCCGCCACATGGAAGCTG", kmer_dict, protein_dict)
    assert(("Q6GZX4", 0) not in occ)

    print('All passed')


if __name__ == "__main__":
    main()
