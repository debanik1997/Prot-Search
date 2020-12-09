import pickle
import numpy as np

def main():
  for i in [1,2,3,4,5,10]:
    f = 'kmer_dict_k_' + str(i) + '_num_prots_100.pickle'
    with open(f, "rb") as pf: kmers = pickle.load(pf)
    print(f'Length l = {i}\nNumber of keys is {len(kmers.keys())}')
    avg = np.mean([len(kmers[k]) for k in kmers.keys()])
    print(f'Average vals per key is {avg}')

main()
