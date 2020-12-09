'''
Python program to crete kmers from uniprot fasta file

Code used for Computational Genomics Projects Fa 2020

Usage:
python kmers.py --infile ../data/protein_dict_num_prots_100.pickle -k 5 -n 100
'''

import argparse
import pickle

OUTDIR = '../data/'

def read_from_pickle(f):
    """ Reads a dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        dictionary = pickle.load(fd)
    return dictionary

def read_pickle(protein_dict, k, n):
    kmers = {} # dictionary to store kmers

    count_proteins = 0 # count of proteins read
    for protein_id in protein_dict.keys():

        # get protein sequence
        protein_seq = protein_dict[protein_id]

        # generate all possible kmers
        for i in range(len(protein_seq) - k + 1):
            kmer = protein_seq[i:i+k] # could be made faster by windowing
            if kmer not in kmers:
                kmers[kmer] = [(protein_id, i)]
            else:
                kmers[kmer].append((protein_id, i))
        count_proteins +=1
        if n == count_proteins:
            break
    return kmers

def read_fasta(lines, k, n, delim):
    print('Making dict ..')

    seq = ''
    kmers = {} # dictionary to store kmers

    count_proteins = 0 # count of proteins read
    for j, line in enumerate(lines):
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
            if n == count_proteins:
                break

            # store next protein id
            protein_id = line.split(delim)[1]

            # initialize seq again
            seq = ''
        else:
            # add to current protein sequence
            seq += line
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
        extension = infile.split('.')[-1]


        if  extension == 'pickle':
            database = read_from_pickle(infile)
        elif extension == 'fasta':
            with open(infile) as f:
                database = f.readlines()
        else:
            print('Please give a .pickle or .fasta infile')
            exit()

        global OUTDIR

        if n == None or n>len(database):
            n = len(database)
            outfile = OUTDIR + 'kmer_dict_k_' + str(k) + '.pickle'
        else:
            outfile = OUTDIR + 'kmer_dict_k_' + str(k) + '_num_prots_' + str(n) + '.pickle'

        if extension == 'fasta':
            kmers = read_fasta(database, k, n, delim)
        else:
            kmers = read_pickle(database, k, n)

        # now store kmers to file using pickle
        print('Dumping dictionary ..')
        # print(len(kmers.keys()))
        with open(outfile, 'wb') as f:
            # HIGHEST_PROTOCOL means use the fastest protocal available
            pickle.dump(kmers, f, pickle.DEFAULT_PROTOCOL)


main()
