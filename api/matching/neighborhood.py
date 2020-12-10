'''
Code to make neighbors

For each offset in protein sequence, try substituting all possible amino acids
and calculating scores

Usage example:
python3 neighborhood.py -i TTSNT -m 1 -t 17

for example: for protein AB, and possible amino acids A,B,F do
'''

import pickle
import argparse
from collections import defaultdict

BLOSUM_FILE = '../../data/BLOSUM62.pickle'

def read_from_pickle(f):
    """ Reads a dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        dictionary = pickle.load(fd)
    return dictionary



def main():

    ap = argparse.ArgumentParser()
    ap.add_argument('--input', '-i', type = str, required = True,
                    help = 'original string to find neighbors for')
    ap.add_argument('-l', type = int, default = 4,
                    help = 'length of lmer')
    ap.add_argument('-m', '--max_mismatches', default = 1, type = int,
                    help = 'maximum mismatches allowed')
    ap.add_argument('-t', default = 11, type = int,
                    help = 'score threshold for keeping neighbors')

    args = ap.parse_args()

    max_mismatch = args.max_mismatches
    threshold = args.t
    l = args.l #kmer length
    input_seq = args.input.strip()

    blosum62 = read_from_pickle(BLOSUM_FILE)
    # make dict from 'columns' key
    blosum_indices = {c:i for i,c in enumerate(blosum62['columns'])}
    possible_acids = set(blosum_indices.keys())

    print('Neighbors\tScore\tMismatch\n')
    # for each lmer in input, introduce max_mismatch errors
    seeds = defaultdict(list)
    # for each lmer in input, introduce max_mismatch errors
    for i in range(len(input_seq) - l + 1):

        # offset for each partition
        part_off = i
        # lmer
        part = input_seq[part_off:part_off+l]

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

    for seed in seeds:
        print(f'{seed}\t\t{seeds[seed][0][2]}\t{seeds[seed][0][1]}')


main()
