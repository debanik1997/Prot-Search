import sys
from six_frame_translation import six_frame_translation
from collections import defaultdict
from utils import *

KMER_FILE = 'kmer_dict_k_4_num_prots_2.pickle'
UNIPROT_FILE = 'uniprot_sprot.fasta'

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

    # Result stored as a set of tuples (id, offset)
    occurrences = set()

    for protein in prot_translations:
        for i in range(60//4): # for all 15 partitions

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
    return occurrences

# Test Cases
def main():

    p = "ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG"

    # 1 Mismatch - should pass
    occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG")
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
