'''
Script to create a map from {protein_id: protein_sequence}

Example:
python3 proteins.py --infile uniprot_sprot.fasta --num_proteins 100
'''
from utils import fasta_iterator
import argparse
import pickle

def fasta_iterator(f, num_proteins, delim='|'):
    """ An iterator to iterate through a fasta file. Yields a tuple of (protein id, protein sequence) """
    with open(f) as fd:
      lines = fd.readlines()
    line_idx = 0

    for i in range(num_proteins):
        line = lines[line_idx].strip()
        protein_id = line.split(delim)[1]
        seq = ""
        line_idx += 1

        while True:
            line = lines[line_idx].strip()
            if line[0] < 'A' or line[0] > 'Z':
                break
            seq += line
            line_idx += 1
        
        yield (protein_id, seq)

def create_protein_map(f, num_proteins):
    protein_dict = {}
    iterator = fasta_iterator(f, num_proteins)
    for (protein_id, seq) in iterator:
        assert(protein_id not in protein_dict)
        protein_dict[protein_id] = seq
    assert(len(protein_dict.keys()) == num_proteins)
    return protein_dict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--infile', type=str, required=True,
                    help='path to FASTA file')
    ap.add_argument('--num_proteins', type=int, default=100,
                    help='Number of proteins to include in the map')

    args = ap.parse_args()
    infile = args.infile
    num_proteins = args.num_proteins
    protein_dict = create_protein_map(infile, num_proteins)

    outfile = f"protein_dict_num_prots_{num_proteins}.pickle"
    with open(outfile, 'wb') as f:
        pickle.dump(protein_dict, f, pickle.DEFAULT_PROTOCOL)


if __name__ == "__main__":
    main()