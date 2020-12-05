import argparse
from index_assisted_dp import query
from pigeonhole_principle_approximate_match import query
from utils import read_from_pickle

def run_test_suite(infile, kmer_dict, protein_dict):
    pass

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--infile', type=str, required=True,
                    help='path to test_reads file')
    ap.add_argument('--kmer_dict_file',type=str, required=True,
                    help='path to kmer dictionary pickle file')
    ap.add_argument('--protein_dict_file', type=std, required=True,
                    help='path to protein dictionary pickle file')
    args = ap.parse_args()
    infile = args.infile
    kmer_dict = read_from_pickle(args.kmer_dict_file)
    protein_dict = read_from_pickle(args.protein_dict_file)

    
    run_test_suite(infile, kmer_dict, protein_dict)

if __name__ == "__main__":
    main()