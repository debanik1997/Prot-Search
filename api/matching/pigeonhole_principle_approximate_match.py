import sys
from six_frame_translation import six_frame_translation
from collections import defaultdict
from utils import *
import pickle
from memory_profiler import profile

KMER_DICT_FILE = 'kmer_dict_k_4_num_prots_100.pickle'
PROTEIN_DICT_FILE = 'protein_dict_num_prots_100.pickle'

def query(p, kmer_dict, protein_dict, max_mismatch=4):
    """ Function that outputs (approximate) matches of a DNA read against 
    a protein database. Implemented using k-mers and the pigeonhole principle.

    Does not account for gaps or insertions (only mismatches)
    """
    min_rl = 4*(max_mismatch + 1)*3
    assert(len(p) >= min_rl, f"Enter a read of length {min_rl} or higher")
    translations = six_frame_translation(p)
    
    # Result stored as a dictionary of (id, offset) : num_mismatches
    occurrences = {}

    for translation in translations:
        for i in range(len(translation) // 4):
            off = 4*i
            part = translation[off:off+4]

            if part not in kmer_dict:
                continue
            
            hits = kmer_dict[part]
            # hits is a list of [(protein_id, offset)]
            for (id, offset) in hits:
                t = protein_dict[id]
                assert(len(t) != 0)
                nmm = 0
                if offset - off + len(translation) <= len(t) and offset - off >= 0:
                    for j in range(0, len(translation)): 
                        if t[offset-off+j] != translation[j]:
                            nmm += 1
                            if nmm > max_mismatch:
                                break
                    if nmm <= max_mismatch:
                        occurrences[(id, offset-off)] = nmm
    return occurrences

def main():
    kmer_dict = read_from_pickle(KMER_DICT_FILE)
    protein_dict = read_from_pickle(PROTEIN_DICT_FILE)

    p = "ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG"

    # 1 Mismatch - should pass
    occ = query("ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCT", kmer_dict, protein_dict)
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

    print("All passed")


if __name__ == "__main__":
    main()