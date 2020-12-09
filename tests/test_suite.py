"""
Prints a summary of execution time and memory usage

Example
python3 test_suite.py --infile test_reads.txt
"""

import argparse
import tracemalloc
import pickle
from timeit import default_timer as timer
import sys
sys.path.append('../api/matching')
from index_assisted_dp import query as index_assisted_query
from pigeonhole_principle_approximate_match import query as pigeonhole_query

KMER_FILE = '../data/kmer_dict_k_4_num_prots_100.pickle'
PROTEIN_FILE = '../data/protein_dict_num_prots_100.pickle'

def read_from_pickle(f):
    """ Reads a dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        dictionary = pickle.load(fd)
    return dictionary

def profile(func):
    def wrapper(reads, kmer_dict, protein_dict):
        tracemalloc.start()
        start_time = timer()
        func(reads, kmer_dict, protein_dict)
        end_time = timer()
        _, peak = tracemalloc.get_traced_memory()
        print(f"Peak Memory Usage: {peak / 10**3}KB")
        avg_time = (end_time - start_time) / len(reads) * 10**3
        print(f"Average Execution Time: {avg_time:.2f}ms \n")
    
    return wrapper

def print_occurences(read, occ):
    ordered = {k: v for k, v in sorted(occ.items(), key=lambda item: item[1])}
    for k, v in ordered.items():
        protein_id, offset = k
        print(f"{protein_id} \t {offset} \t \t {v}")

@profile
def pigeonhole_test(reads, kmer_dict, protein_dict):
    print("l-mer Filtration Algorithm:")
    for read in reads:
        occ = pigeonhole_query(read, kmer_dict, protein_dict)
        if len(occ) == 0:
            print("No occurrences found")
            continue
        print(f"Matches for read {read}")
        print(f"Id \t Offset \t Mismatches")
        print_occurences(read, occ)

@profile
def index_assisted_dp_test(reads, kmer_dict, protein_dict):
    print("Index Assisted DP Algorithm:")
    for read in reads:
        occ = index_assisted_query(read, kmer_dict, protein_dict)
        if len(occ) == 0:
            print("No occurrences found")
            continue
        print(f"Matches for read {read}")
        print(f"Id \t Offset \t Mismatches/Gaps/Insertions")
        print_occurences(read, occ)

def run_test_suite(infile, kmer_dict, protein_dict):
    with open(infile) as f:
        reads = f.readlines()
    reads = [read.strip("\n") for read in reads]
    pigeonhole_test(reads, kmer_dict, protein_dict)
    index_assisted_dp_test(reads, kmer_dict, protein_dict)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--infile', type=str, required=True,
                    help='path to test_reads file')
    args = ap.parse_args()
    infile = args.infile
    kmer_dict = read_from_pickle(KMER_FILE)
    protein_dict = read_from_pickle(PROTEIN_FILE)

    run_test_suite(infile, kmer_dict, protein_dict)

if __name__ == "__main__":
    main()