'''
Prints summary of creation runtime and size of dictionary for l-mer preprocessing

Usage:
python3 kmer_dict_tests.py
'''

import pickle
import numpy as np
import os
import time

def main():
  DATA_PATH = '../data/'
  for i in [1,2,3,4,5,10]:
    f = DATA_PATH + 'kmer_dict_k_' + str(i) + '_num_prots_100.pickle'
    with open(f, "rb") as pf: kmers = pickle.load(pf)
    print(f'Length l = {i}\nNumber of keys is {len(kmers.keys())}')
    avg = np.mean([len(kmers[k]) for k in kmers.keys()])
    print(f'Average vals per key is {avg}\n')


  print('\n\nTime taken to make each lmer dictionary, with number of proteins = 100')

  for i in [1,2,3,4,5,10]:

      cmd = f'python ../scripts/kmers.py --infile ../data/protein_dict_num_prots_100.pickle -k {i} -n 100 > output'
      time.sleep(1)
      start_time = time.time()
      for _ in range(2): # run only 2 times -> higher than this makes pickle fail
          # run this command
          os.system(cmd)
      end_time = time.time()
      print(f'Time taken to make dictionary for l = {i}:\t {(end_time - start_time)/2} s')


  print('\n\nTime taken to make dictionary for different number of proteins, l = 4')

  for i in [10,100,200,500,1000]:

      cmd = f'python ../scripts/kmers.py --infile ../data/protein_dict_num_prots_1000.pickle -k 4 -n {i} > output'
      time.sleep(1)
      start_time = time.time()
      for _ in range(10): # run 10 times
            # run this command
            os.system(cmd)
      end_time = time.time()
      print(f'Time taken to make dictionary for n = {i}:\t {(end_time - start_time)/10} s')
      size = os.path.getsize(DATA_PATH + 'kmer_dict_k_4_num_prots_' + str(i) + '.pickle')
      print(f'\t\t\tSize of file:\t{size/1024} KB')


main()
