"""
Script to test and benchmark reads
Prints a summary of execution time and memory usage

Example
python3 test_suite.py --infile test_reads.txt --kmer_dict_file kmer_dict_k_4_num_prots_100.pickle --protein_dict_file protein_dict_num_prots_100.pickle 
"""

import argparse
import tracemalloc
from timeit import default_timer as timer
from index_assisted_dp import query as index_assisted_query
from pigeonhole_principle_approximate_match import query as pigeonhole_query
from utils import read_from_pickle

def print_occurences(read, occ):
    print(f"Matches for read {read}")
    print(f"Id \t Offset \t Mismatches/Gaps/Insertions")
    ordered = {k: v for k, v in sorted(occ.items(), key=lambda item: item[1])}
    for k, v in ordered.items():
        protein_id, offset = k
        print(f"{protein_id} \t {offset} \t \t {v}")

def pigeonhole_test(reads, kmer_dict, protein_dict):
    print("l-mer Filtration Technique:")
    tracemalloc.start()
    start_time = timer()
    for read in reads:
        occ = pigeonhole_query(read, kmer_dict, protein_dict)
        print_occurences(read, occ)
    end_time = timer()
    _, peak = tracemalloc.get_traced_memory()
    print(f"Peak Memory usage by l-mer filtration algorithm is {peak / 10**3}KB")
    print(f"Average Execution time by l-mer filtration algorithm is {(end_time-start_time)/len(reads) * 10**3}ms")
    print()

def index_assisted_dp_test(reads, kmer_dict, protein_dict):
    tracemalloc.start()
    start_time = timer()
    for read in reads:
        occ = index_assisted_query(read, kmer_dict, protein_dict)
        print_occurences(read, occ)
    end_time = timer()
    _, peak = tracemalloc.get_traced_memory()
    print(f"Peak Memory usage by index assisted query algorithm is {peak / 10**3}KB")
    print(f"Average Execution time by index assisted query algorithm is {(end_time-start_time)/len(reads) * 10**3}ms")
    print()


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
    ap.add_argument('--kmer_dict_file',type=str, required=True,
                    help='path to kmer dictionary pickle file')
    ap.add_argument('--protein_dict_file', type=str, required=True,
                    help='path to protein dictionary pickle file')
    args = ap.parse_args()
    infile = args.infile
    kmer_dict = read_from_pickle(args.kmer_dict_file)
    protein_dict = read_from_pickle(args.protein_dict_file)

    run_test_suite(infile, kmer_dict, protein_dict)

if __name__ == "__main__":
    main()