'''
Python program to create kmers from uniprot fasta file

Code used for Computational Genomics Projects Fa 2020

Usage:
python kmers.py --infile uniprot_sprot.fasta -k 5 -n 100
'''

import argparse
import pickle

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
    with open(infile) as f:
      lines = f.readlines()

    if n == None:
      n = len(lines)
      outfile = 'kmer_dict_k_' + str(k) + '.pickle'
    else:
      outfile = 'kmer_dict_k_' + str(k) + '_num_prots_' + str(n) + '.pickle'


    print('Making dict ..')
    kmers = {} # dictionary to store kmers

    count_proteins = 0 # count of proteins read
    for line in lines:
      line = line.strip()
      if line[0] <'A' or line[0] > 'Z': # id line of protein
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

    # now store kmers to file using pickle
    print('Dumping dictionary ..')
    with open(outfile, 'wb') as f:
      # HIGHEST_PROTOCOL means use the fastest protocal available
      pickle.dump(kmers, f, pickle.HIGHEST_PROTOCOL)


main()
