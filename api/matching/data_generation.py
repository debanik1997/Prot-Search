'''
Script to generate random reads from a given FASTA file

Example Usage:
python3 data_generation.py --infile uniprot_sprot.fasta --seed 1 --num_errors 4 --num_reads 5 --read_length 60 

Writes reads to an output file test_reads.txt
'''
import argparse
import random
from utils import fasta_iterator

# Duplicate keys would be overwritten with a value that came later
inverse_codon_map = {
    'F':'UUU',
    'F':'UUC',
    'L':'UUA',
    'L':'UUG',
    'S':'UCU',
    'S':'UCC',
    'S':'UCA',
    'S':'UCG',
    'Y':'UAU',
    'Y':'UAC',
    '*':'UAA',
    '*':'UAG',
    'C':'UGU',
    'C':'UGC',
    '*':'UGA',
    'W':'UGG',
    'L':'CUU',
    'L':'CUC',
    'L':'CUA',
    'L':'CUG',
    'P':'CCU',
    'P':'CCC',
    'P':'CCA',
    'P':'CCG',
    'H':'CAU',
    'H':'CAC',
    'Q':'CAA',
    'Q':'CAG',
    'R':'CGU',
    'R':'CGC',
    'R':'CGA',
    'R':'CGG',
    'I':'AUU',
    'I':'AUC',
    'I':'AUA',
    'M':'AUG',
    'T':'ACU',
    'T':'ACC',
    'T':'ACA',
    'T':'ACG',
    'N':'AAU',
    'N':'AAC',
    'K':'AAA',
    'K':'AAG',
    'S':'AGU',
    'S':'AGC',
    'R':'AGA',
    'R':'AGG',
    'V':'GUU',
    'V':'GUC',
    'V':'GUA',
    'V':'GUG',
    'A':'GCU',
    'A':'GCC',
    'A':'GCA',
    'A':'GCG',
    'D':'GAU',
    'D':'GAC',
    'E':'GAA',
    'E':'GAG',
    'G':'GGU',
    'G':'GGC',
    'G':'GGA',
    'G':'GGG', 
}

codon_map = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
}

def generate_test_reads(infile, seed, num_errors, num_reads, n):
    # Initialize random seed
    random.seed(seed)

    # We want a protein sequence of n/3 since 3 DNA Codons map to 1 protein
    p = n // 3

    # Outfile
    f = open('test_reads.txt', 'w+')

    # Iterator
    iterator = fasta_iterator(infile, num_reads)

    for (protein_id, seq) in iterator:
        # Find a random read of length p
        start_idx = random.randint(0, len(seq)-p-1)
        end_idx = start_idx + p
        rand_prot = seq[start_idx:end_idx]

        # Expand the rand_seq to a DNA Read using the inverse codon map
        rand_dna = []
        for c in rand_prot:
            dna = inverse_codon_map[c].replace("U", "T")
            rand_dna += list(dna)

        # TODO: Insert num_errors errors in rand_dna. Currently only mismatches
        for i in range(num_errors):
            error_idx = random.randint(0, len(rand_dna)-1)
            # Replace with a random codon
            old_codon = rand_dna[error_idx]
            rand_idx = random.randint(0, 3)
            rand_dna[error_idx] = codon_map[rand_idx]

        rand_dna = "".join(rand_dna)
        # Write to out-file
        f.write(f"{rand_dna}\n")

    f.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--infile', type = str, required = True,
                        help = 'path to fastA file')
    ap.add_argument('--seed', type=int, default=1, 
                    help='SEED value to generate reads')
    ap.add_argument('--num_errors', type=int, default=0,
                    help='Maximum number of errors (mismatches/gaps/insertions) to insert')
    ap.add_argument('--num_reads', type=int, default=5,
                    help='Maximum number of reads to generate')
    ap.add_argument('--read_length', type=int, default=60,
                    help='Read length')

    args = ap.parse_args()
    infile = args.infile
    seed = args.seed
    num_errors = args.num_errors
    num_reads = args.num_reads
    n = args.read_length

    generate_test_reads(infile, seed, num_errors, num_reads, n)

if __name__ == "__main__":
    main()