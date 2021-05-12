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


## Needed Libraries

Biopython: https://biopython.org/
