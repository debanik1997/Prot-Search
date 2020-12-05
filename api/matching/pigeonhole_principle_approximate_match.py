import sys
from six_frame_translation import six_frame_translation
from collections import defaultdict
from utils import *
import pickle

KMER_FILE = 'kmer_dict_k_4_num_prots_2.pickle'
UNIPROT_FILE = 'uniprot_sprot.fasta'

def query(p, kmer_dict, protein_dict, max_mismatch=4):
    """ Function that outputs (approximate) matches of a DNA read against 
    a protein database. Implemented using k-mers and the pigeonhole principle.

    Default number of mismatches is fixed (4)
    Assumes length of p is fixed (60)
    Does not account for gaps or insertions (only mismatches)
    """
    assert (len(p) == 60)

    translations = six_frame_translation(p)
    
    # Result stored as a set of tuples (id, offset)
    occurrences = set()

    for translation in translations:
        for i in range(5):
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
                        occurrences.add((id, offset-off))  
    return occurrences

# Test Cases
def main():
     # kmer_dict holds mapping from {kmer: [(protein_id, offset)]}
    kmer_dict = read_from_pickle(KMER_FILE)

    # protein_dict holds mapping from {protein_id : protein_seq}
    protein_dict = create_protein_map(UNIPROT_FILE, 1)# Store data (serialize)

    p = "ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG"

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


if __name__ == "__main__":
    main()