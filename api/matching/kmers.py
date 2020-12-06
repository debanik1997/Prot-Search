'''
Python program to crete kmers from uniprot fasta file

Code used for Computational Genomics Projects Fa 2020

Usage:
python kmers.py --infile protein_dict_num_prots_100.pickle -k 4 -n 100
'''

import argparse
import pickle
from utils import read_from_pickle

# to be used with uniprot pickle files
PROTEIN_DICT_FILE = 'protein_dict_num_prots_100.pickle'

def dict_from_fasta(lines:list, max_prots:int, k:int) ->dict:
    print('Making dict ..')
    kmers = {}
    count_proteins = 0 # count of proteins read
    for i, line in enumerate(lines):
        line = line.strip()
        if (line[0] <'A' or line[0] > 'Z'): # id line of protein
            # if you've hit the next protein seq, that means you have the
            # complete last one. Make kmers for it
            for i in range(len(seq) - k + 1):

                kmer = seq[i:i+k] # could be made faster by windowing
                if kmer not in kmers:
                    kmers[kmer] = [(protein_id, i)]
                else:
                    kmers[kmer].append((protein_id, i))

            count_proteins +=1
            if max_prots == count_proteins:
                break

            # store next protein id
            protein_id = line.split(delim)[1]

            # initialize seq again
            seq = ''
        else:
            # add to current protein sequence
            seq += line
    return kmers

def dict_from_pickle(proteins:dict, max_prots:int, k:int)->dict:

    print('Making dict ..')
    kmers = {}

    count_proteins = 0 # count of proteins read
    for protein_id in proteins.keys():
        seq = proteins[protein_id]

        for i in range(len(seq) - k + 1):

            kmer = seq[i:i+k] # could be made faster by windowing
            if kmer not in kmers:
                kmers[kmer] = [(protein_id, i)]
            else:
                kmers[kmer].append((protein_id, i))

        count_proteins +=1
        if max_prots == count_proteins:
            break

    return kmers

def main():
        ap = argparse.ArgumentParser()
        ap.add_argument('--infile', type = str, required = True,
                        help = 'path to fastA file')
        ap.add_argument('-k', type = int, default = 4,
                        help = 'length of kmer')
        ap.add_argument('-d', '--delimiter', default = '|',
                        help = 'delimiter in fasta file')
        ap.add_argument('-n', default = None, type = int,
                        help = 'maximum number or proteins to read from file')

        args = ap.parse_args()

        k = args.k
        infile = args.infile

        delim = args.delimiter
        n = args.n

        seq = ''

        # try loading the fasta file
        if infile.split('.')[-1] == 'fasta':
            with open(infile) as f:
                proteins = f.readlines()
            kmers = dict_from_fasta(proteins, n, k)
        elif infile.split('.')[-1] == 'pickle':
            # load pickle file
            # protein_dict holds mapping from {protein_id : protein_seq}
            proteins = read_from_pickle(PROTEIN_DICT_FILE)
            kmers = dict_from_pickle(proteins, n, k)
        else:
            print(f"Invalid file extension {infile.split('.')[-1]}. Try .fasta or .pickle")
        if n == None:
            n = len(proteins) # works for list and dict
            outfile = 'kmer_dict_k_' + str(k) + '.pickle'
        else:
            n = min(n, len(proteins)) # in case n is less than actual # of proteins
            outfile = 'kmer_dict_k_' + str(k) + '_num_prots_' + str(n) + '.pickle'

        # now store kmers to file using pickle
        print('Dumping dictionary ..')
        # print(kmers)
        with open(outfile, 'wb') as f:
            # HIGHEST_PROTOCOL means use the fastest protocal available
            pickle.dump(kmers, f, pickle.DEFAULT_PROTOCOL)


main()
