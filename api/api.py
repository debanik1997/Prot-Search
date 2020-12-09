#----------------------------------------------------------------------------#
# Imports
#----------------------------------------------------------------------------#
import sys
sys.path.append('./matching')
import pigeonhole_principle_approximate_match
import index_assisted_dp
from flask import Flask, request
import time
import os
import pickle

KMER_FILE = '../data/kmer_dict_k_4_num_prots_100.pickle'
UNIPROT_FILE = '../data/protein_dict_num_prots_100.pickle'

#----------------------------------------------------------------------------#
# Functions
#----------------------------------------------------------------------------#

def read_from_pickle(f):
    """ Reads a dictionary from file """
    kmer_dict = {}
    with open(f, "rb") as fd:
        dictionary = pickle.load(fd)
    return dictionary

#----------------------------------------------------------------------------#
# App Config.
#----------------------------------------------------------------------------#

app = Flask(__name__)

#----------------------------------------------------------------------------#
# Controllers.
#----------------------------------------------------------------------------#


@app.route('/api/time')
def get_current_time():
    return {'time': time.time()}


@app.route('/api/protein', methods=['GET', 'POST'])
def protein_query():
    if request.method == 'POST':
        data = request.get_json()
        kmer_dict = read_from_pickle(KMER_FILE)
        protein_dict = read_from_pickle(UNIPROT_FILE)
        if data['gaps_allowed']:
            return {'match': list( index_assisted_dp.query(data['pattern'], kmer_dict, protein_dict, data["max_mismatches"]))}
        else:
            return {'match': list(pigeonhole_principle_approximate_match.query(data['pattern'], kmer_dict, protein_dict, data["max_mismatches"]))}
    else:
        return {"status": "ok"}

#----------------------------------------------------------------------------#
# Launch.
#----------------------------------------------------------------------------#


# Default port:
if __name__ == '__main__':
    app.run()
