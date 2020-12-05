#----------------------------------------------------------------------------#
# Imports
#----------------------------------------------------------------------------#

import sys
from flask import Flask, request
import time
import os
from matching.pigeonhole_principle_approximate_match import query

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
        pattern = data["pattern"]
        return {"pattern": pattern}
    else:
        return {"Method": "Get"}

#----------------------------------------------------------------------------#
# Launch.
#----------------------------------------------------------------------------#

# Default port:
if __name__ == '__main__':
    app.run()
