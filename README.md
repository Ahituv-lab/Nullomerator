===========
# Nullomerator
===========

## Installation

1.- Clone github repository

	git clone https://github.com/Ahituv-lab/Nullomerator

## Functions

**1.EnumerateNullomers**
Extracts all nullomers of specified kmer lengths in a FASTA sample.

**2.ExtractMutationNullomers**
Finds all mutations that cause the resurfacing of a list of nullomers.

**3.IdentifyRecurrentNullomers**
Identifies nullomers that recur in a dataset through mutagenesis.

**4.FindAlmostNullomers**
Identifies the positions that can create a list of nullomers genome-wide, for every possible substitution and single base-pair insertion and deletion.

**5.FindNullomerVariants**
Removes nullomers that are likely to result from common variants in a user specified variant VCF file

**6.FindDNANullomersFromReads**
Performs the identification of nullomers in raw read samples


## Third Party Resources

Biopython: https://biopython.org/

## Documentation (Draft, to assist with testing)

### EnumerateNullomers

The command for running EnumerateNullomers is as follows:
```
python3 EnumerateNullomers.py --genome_file --nullomer_output_filepath --kmer_output_filepath --kmer_length
```

#### `--genome_file`
Path to the fasta file containing the genome being analyzed. This file must be formatted such that the FASTA headers read ">chr1", ">chr2", "chr3", ... , "chrX", "chrY".

#### `--nullomer_output_filepath`
Path to the output .txt file where the nullomers absent in the supplied genome will be written
