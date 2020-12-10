'''
This python code matches given DNA read with protein database using a kmer index
and a dp based edit distance

Written for Group 63 Project, Computational Genomics Fa 2020

Usage:

python3 blast.py

'''

import sys
from six_frame_translation import six_frame_translation
from collections import defaultdict
from utils import *
import numpy as np

KMER_DICT_FILE = 'kmer_dict_k_4_num_prots_2.pickle'
PROTEIN_DICT_FILE = 'protein_dict_num_prots_100.pickle'
UNIPROT_FILE = 'uniprot_sprot.fasta'
BLOSUM_FILE = 'BLOSUM62.pickle'

blosum62 = {}
blosum_indices = {}
possible_acids = set()

def calculate_neighbors(input, max_mismatch = 4, l = 4, threshold = 10)->list:
    global blosum62, blosum_indices, possible_acids
    if blosum62 is None:
        # load once and store

        blosum62 = read_from_pickle(BLOSUM_FILE)
        # make dict from 'columns' key
        blosum_indices = {c:i for i,c in enumerate(blosum62['columns'])}
        possible_acids = set(blosum_indices.keys())

    seeds = defaultdict(list)
    # for each lmer in input, introduce max_mismatch errors
    for i in range(len(input) - l + 1):

        # offset for each partition
        part_off = i
        # lmer
        part = input[part_off:part_off+l]

        # will store substrings and score for last iteration
        last = [['',0,0]]
        # last and current are lists of format ([str, score, mismatches])

        # will store substrings and score for current iteration
        cur = []

        # for all amino acids we can substitute in this lmer
        for length in range(l):
          for prev in last:
            # for all substrings for last iteration, try concatenating all
            # possible amino acids and adding score from this iteration

            # amino acid in input we are aligning with
            actual_a = part[length]

            # if we have hit max mismatch threshold, copy the rest of protein as is
            if prev[2] == max_mismatch:
                # just append the correct protein, and move on
                cur.append([prev[0] + actual_a,
                            prev[1] + int(blosum62[actual_a][blosum_indices[actual_a]]),
                            prev[2]])
                continue

            # otherwise, try finding other neighbors
            for next_a in possible_acids:
              delta = 1
              if next_a == actual_a:
                delta = 0
              cur.append([prev[0] + next_a,
                        prev[1] + int(blosum62[actual_a][blosum_indices[next_a]]),
                        prev[2] + delta])


          last = cur
          # empty current for next iteration
          cur = []

        # remove all below the threshold
        for prev in last:
            if prev[1]>=threshold:
                # seeds has format: seed -> [[offset, mismatches, score]]
                seeds[prev[0]].append([part_off, prev[2], prev[1]])
    return seeds


def query(read, kmer_dict, protein_dict, max_mismatch = 4, l = 4):
    """ Function that outputs (approximate) matches of a DNA read against
    a protein database. Implemented using k-mers and the pigeonhole principle.

    Default value of l-mer length is 4
    Default value of max mismatches allowed is 4

    """

    # will throw AssertionError if we can't make >max_mismatch partitions of size l
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

    score_threshold = 5
    hsp = []

    for protein in prot_translations:
        seeds = calculate_neighbors(protein)

        for seed in seeds.keys(): # for all seeds
            if seed not in kmer_dict:
                continue

            # hits is a list of [(protein_id, offset_in_protein)]
            hits = kmer_dict[seed]

            for (id, offset) in hits:

                t = protein_dict[id]

                if len(t) == 0:
                    continue


                # seeds has format: seed -> [[offset, mismatches, score]]
                # seeds[seed] is a list of all the places 'seed' occured in read

                for seed_parts in seeds[seed]:
                    # offset for each partition in read
                    part_off = seed_parts[0]


                    # start from here, go left and right.
                    # for this implementation, choose what to do first based
                    # on scores on left and right

                    # if |max_score - cur_score| falls below a threshold, stop
                    max_score = seed_parts[2]
                    left_score, right_score = max_score, max_score

                    # initialize left and right with offset of partition edges
                    lf_seed, rt_seed = part_off, part_off + l -1

                    # initialize left and right for database protein
                    lf_t, rt_t = off, off + l - 1
                    # bools for expanding left and right
                    go_left, go_right = True, True
                    while go_left and go_right:

                        if go_left:
                            lf_seed, lf_t = lf_seed - 1, lf_t - 1
                            if left < 0: # cannot go left anymore
                                go_left = False
                                lf_score = 0
                                lf_t +=1

                            else:
                                actual_a, read_a = t[lf_t], protein[lf_seed]

                                # add score
                                lf_score = int(blosum62[actual_a][blosum_indices[cur_a]])

                        if go_right:
                            rt_seed, rt_t = rt_seed + 1, rt_t + 1
                            if rt_seed >= len(protein) or rt_t >= len(t):
                                go_right = False
                                rt_score = 0
                                rt_t -=1
                            else:
                                actual_a, read_a = t[rt_t], protein[rt_seed]
                                rt_score = int(blosum62[actual_a][blosum_indices[cur_a]])

                        running_score += lf_score, rt_score # extend in both directions
                        max_score = max(running_score, max_score)
                        diff_from_max = max_score - running_score

                        if diff_from_max >= score_threshold:
                            continue
                        else:
                            # below threshold: abandon
                            # check if extending in only one direction is better

                            if rt_score > lf_score and  diff_from_max + lf_score >= score_threshold:
                                # if left score was smaller, and removing it
                                # brings total score above thresh, stop extending left
                                go_left = False
                                lf_seed, lf_t = lf_seed -1, lf_t -1
                                running_score -= lf_score
                                lf_score = 0
                            elif diff_from_max + rt_score >= score_threshold:
                                # stop extending right
                                go_right = False
                                rt_seed, rt_t = rt_seed -1, rt_t-1
                                running_score -= rt_score
                                rt_score = 0
                            else:
                                # stop complete alignment
                                go_left = False
                                go_right = False

                    # out of while loop
                    if rt_t - lf_t > l: # if there was some extension
                        hsp.append([t, lf_t, rt_t])


                # go through all HSP's now
                # todo
                if (id, lf, rt) in seen:
                    continue

                seen.add((id, lf, rt))
                raise NotImplementedError

    return occurrences

# Test Cases
def main():

    # kmer_dict holds mapping from {kmer: [(protein_id, offset)]}
    kmer_dict = read_from_pickle(KMER_DICT_FILE)

    # protein_dict holds mapping from {protein_id : protein_seq}
    protein_dict = read_from_pickle(PROTEIN_DICT_FILE)

    global blosum62, blosum_indices, possible_acids
    # loading blosum62 dict
    blosum62 = read_from_pickle(BLOSUM_FILE)
    # make dict from 'columns' key
    blosum_indices = {c:i for i,c in enumerate(blosum62['columns'])}

    possible_acids = set(blosum_indices.keys())

    #########################################
    # Testing on some inputs
    #########################################
    # test protein = "ATGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG"
    # pid = "Q6GZX4"

    # length = 50 - should fail
    try:
        occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCG", kmer_dict, protein_dict, blosum62)
    except AssertionError as e:
        # will print assertion error
        print(e, '. Continuing with other tests.')

    # 1 insertion and 1 deletion - should pass
    occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGCG", kmer_dict, protein_dict, blosum62)
    assert ("Q6GZX4", 0) in occ, print(occ)

    # 1 Mismatch - should pass
    occ = query("TTGGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict, blosum62)
    assert(("Q6GZX4", 0) in occ)

    # 2 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTGCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict, blosum62)
    assert(("Q6GZX4", 0) in occ)

    # 3 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTTCTGAAAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict, blosum62)
    assert(("Q6GZX4", 0) in occ)

    # 4 Mismatches - should pass
    occ = query("TTTGCGTTTAGCGCGGAAGATGTTCTGAGAGAATATGATCGCCGCCGCCGCATGGAAGCG", kmer_dict, protein_dict, blosum62)
    assert(("Q6GZX4", 0) in occ)

    # > 4 Mismatches - should fail
    occ = query("TTGCTTTAGCGCGGAAGATGTTCTGAGAGAAGTATGATCGCCGCCGCCACATGGAAGCTG", kmer_dict, protein_dict, blosum62)
    assert(("Q6GZX4", 0) not in occ)

    print('All passed')


if __name__ == "__main__":
    main()
