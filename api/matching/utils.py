import pickle

def read_from_pickle(f):
    """ Reads the kmer dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        kmer_dict = pickle.load(fd)
    return kmer_dict

def create_protein_map(f, num_proteins, delim='|'):
    """ Creates a map from {protein_id: protein_sequence} for num_proteins proteins """
    with open(f) as fd:
      lines = fd.readlines()
    
    protein_dict = {}
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
        protein_dict[protein_id] = seq

    return protein_dict