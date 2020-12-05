import pickle

def read_from_pickle(f):
    """ Reads a dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        dictionary = pickle.load(fd)
    return dictionary

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