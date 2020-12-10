# Prot-Search
Computational Genomics Final Project - Shreyas Aiyar, Debanik Purkayastha, Ishita Tripathi


  * [Directory Structure](#directory-structure)
  * [Requirements:](#requirements-)
  * [Web App:](#web-app-)
    + [To start server:](#to-start-server-)
    + [To start client:](#to-start-client-)
  * [Testing](#testing)
    + [Generating Synthetic Reads](#generating-synthetic-reads)
    + [Benchmarking/Profiling Algorithms](#benchmarking-profiling-algorithms)
    + [Testing L-Mer Dictionary](#testing-l-mer-dictionary)
  * [Misc](#misc)
    + [Generating L-mer Dictionary File](#generating-l-mer-dictionary-file)
    + [Generating Protein Dictionary File](#generating-protein-dictionary-file)

## Directory Structure
    .
    ├── api
        ├── matching
            ├── index_assisted_dp.py                        # Index Assisted DP algorithm implementation
            ├── neighbourhood.py                            # Neighbourhood matching algorithm implementation
            ├── pigeonhole_principle_approximate_match.py   # L-mer filtration algorithm implementation
            ├── six_frame_translation.py                    # Six frame translation implementation
        ├── api.py
    ├── client                                              # Web app files
    ├── tests
        ├── data_generation.py                              # Script to generate test reads
        ├── kmer_dict_tests.py                              # Test script to test different k values
        ├── test_suite.py                                   # Test script to run and benchmark algorithms
    ├── scripts
        ├── kmers.py                                        # Script to generate L-mer dictionary file
        ├── proteins.py                                     # Script to generate protein id dictionary file
    ├── data                                                # Stores supporting dictionary files in .pickle format
    ├── contributions.txt
    └── README.md

## Requirements:
```
TODO:
```
## Web App:
### To start server:
From root directory
```
cd api && pip install -r requirements.txt && flask run
```
### To start client:
From root directory
```
cd client && yarn install && yarn start
```
## Testing

### Generating Synthetic Reads
```
cd tests
python3 data_generation.py --infile uniprot_sprot.fasta --seed 1 --num_errors 4 --num_reads 5 --read_length 60 --no-gaps
```

### Benchmarking/Profiling Algorithms

Note: Generate synthetic reads prior to this
```
python3 test_suite.py --infile test_reads.txt --mm 4 --l 4
```

### Testing L-Mer Dictionary
```
TODO:
```

## Misc

### Generating L-mer Dictionary File

```
TODO:
```

### Generating Protein Dictionary File

```
TODO:
```