'''
Script to generate random reads from a given FASTA file

Example Usage (without insertions):
python3 data_generation.py --infile uniprot_sprot.fasta --seed 1 --num_errors 4 --num_reads 5 --read_length 60 --no-gaps

Example Usage (insertions):
python3 data_generation.py --infile uniprot_sprot.fasta --seed 1 --num_errors 4 --num_reads 5 --read_length 60 --gaps

Writes reads to an output file test_reads.txt
'''
import argparse
import random

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

def generate_test_reads(infile, seed, num_errors, num_reads, n, gaps):
    # Initialize random seed
    random.seed(seed)

    keys = list(inverse_codon_map.keys())
    proteins = list(filter(lambda x: x != '*', keys))

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
        rand_prot = list(seq[start_idx:end_idx])

        if gaps:
            choices = ["Mismatch", "Gap", "Insertion"]
        else:
            choices = ["Mismatch"]

        # Insert errors in the protein
        for i in range(num_errors):
            error_idx = random.randint(0, len(rand_prot) - 1)
            # Randomly either replace/insert/delete from this position
            choice = random.choice(choices)
            if choice == "Mismatch":
                # Replace with a random protein
                prot = random.choice(proteins)
                rand_prot[error_idx] = prot
            if choice == "Gap":
                rand_prot.pop(error_idx)
            if choice == "Insertion":
                # Insert a random protein
                prot = random.choice(proteins)
                rand_prot.insert(error_idx, prot)

        # Expand the rand_seq to a DNA Read using the inverse codon map
        rand_dna = []
        for c in rand_prot:
            dna = inverse_codon_map[c].replace("U", "T")
            rand_dna.append(dna)

        rand_dna = "".join(rand_dna)

        if (not gaps):
            assert len(rand_dna) == n

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
    ap.add_argument('--gaps', dest='gaps', action='store_true')
    ap.add_argument('--no-gaps', dest='gaps', action='store_false')
    ap.set_defaults(gaps=False)

    args = ap.parse_args()

    infile = args.infile
    seed = args.seed
    num_errors = args.num_errors
    num_reads = args.num_reads
    n = args.read_length
    gaps = args.gaps

    generate_test_reads(infile, seed, num_errors, num_reads, n, gaps)

if __name__ == "__main__":
    main()