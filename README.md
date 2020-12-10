# Prot-Search
Computational Genomics Final Project - Shreyas Aiyar, Debanik Purkayastha, Ishita Tripathi

## Table of Contents
  * [Table of Contents](#table-of-contents)
  * [Directory Structure](#directory-structure)
  * [Requirements](#requirements)
  * [Web App](#web-app)
  * [Testing](#testing)
  * [Misc](#misc)
## Directory Structure
    .
    ├── api
        ├── matching
            ├── index_assisted_dp.py                        # Index Assisted DP algorithm implementation
            ├── neighbourhood.py                            # Neighbourhood search algorithm implementation
            ├── pigeonhole_principle_approximate_match.py   # L-mer filtration algorithm implementation
            ├── six_frame_translation.py                    # Six frame translation implementation
        ├── api.py
    ├── client                                              # Web app files
    ├── tests
        ├── data_generation.py                              # Script to generate test reads
        ├── kmer_dict_tests.py                              # Test script to compare kmer dictionary sizes and creation runtime for different k's and number of proteins. 
        ├── test_suite.py                                   # Test script to run and benchmark algorithms
    ├── scripts
        ├── kmers.py                                        # Script to generate L-mer dictionary file
        ├── proteins.py                                     # Script to generate protein id dictionary file
    ├── data                                                # Stores supporting dictionary files in .pickle format
    ├── contributions.txt
    ├── recording                                           # Screen recording of the web app
    └── README.md

## Requirements

To run the web app specifically, make sure your machine has Node.js installed. After installing Node.js, make sure that yarn is accessible as a package manager.

Node.js - https://nodejs.org/en/download/  
yarn - https://classic.yarnpkg.com/en/docs/install/#mac-stable  

Also, to generate some supporting files, you would require the Uniprot FASTA file which can be downloaded from:  
Uniprot FASTA File - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz  

Place the Uniprot File in /data directory once downloaded

## Web-App
### To start server
From root directory
```
cd api && pip install -r requirements.txt && flask run
```
### To start client
From root directory
```
cd client && yarn install && yarn start
```
## Testing

### Generating Synthetic Reads
```
cd tests
python3 data_generation.py --infile ../data/uniprot_sprot.fasta --seed 1 --num_errors 4 --num_reads 5 --read_length 60 --no-gaps
cd ..
```

### Benchmarking/Profiling Algorithms

Note: Generate synthetic reads prior to this  
Smaller values of L would take a while to run
```
cd tests
python3 test_suite.py --infile test_reads.txt --mm 4 --l 4
cd ..
```

### Testing L-Mer Dictionary
```
cd tests
python3 kmer_dict_tests.py
cd ..
```

## Misc

### Generating L-mer Dictionary File

```
cd scripts
python kmers.py --infile ../data/protein_dict_num_prots_100.pickle -k 5 -n 100
cd ..
```

### Generating Protein Dictionary File

```
cd scripts
python3 proteins.py --infile ../data/uniprot_sprot.fasta --num_proteins 100
cd ..
```

### Neighborhood search for a sequence [PROT_SEQ]
```
cd api/matching
python neighborhood.py -i [PROT_SEQ] -m 1 -t 17
cd ../../
```
