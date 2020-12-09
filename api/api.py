#----------------------------------------------------------------------------#
# Imports
#----------------------------------------------------------------------------#
import sys
sys.path.append('./matching')

import pickle
import os
import time
from flask import Flask, request, jsonify
import index_assisted_dp
import pigeonhole_principle_approximate_match

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
            matches = index_assisted_dp.query(data['pattern'], kmer_dict, protein_dict, data["max_mismatches"])
            matches_as_list = [(k[0], k[1], v) for k,v in matches.items()] 
            return {"match": matches_as_list}
        else:
            matches = pigeonhole_principle_approximate_match.query(data['pattern'], kmer_dict, protein_dict, data["max_mismatches"])
            matches_as_list = [(k[0], k[1], v) for k,v in matches.items()] 
            return {"match": matches_as_list}
    else:
        return {"status": "ok"}

#----------------------------------------------------------------------------#
# Launch.
#----------------------------------------------------------------------------#


# Default port:
if __name__ == '__main__':
    app.run()
