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
python3 EnumerateNullomers.py --genome_filepath --nullomer_output_filepath --kmer_output_filepath --kmer_length
```
**Main arguments**

###### `--genome_filepath`
Path to the fasta file containing the genome being analyzed. This file must be formatted such that the FASTA headers read ">chr1", ">chr2", "chr3", ... , "chrX", "chrY".

###### `--nullomer_output_filepath`
Path to the output .txt file where the nullomers absent from the supplied genome will be written, with one nullomer written on each line. 

###### `--kmer_output_filepath`
Path to the output .tsv file where the kmers, and their corresponding occurrence counts, from the supplied genome will be written. 

###### `--kmer_length`
Length of kmers/nullomers to be enumerated.


### ExtractMutationNullomers

The command for running ExtractMutationNullomers is as follows:
```
python3 ExtractMutationNullomers.py --genome_filepath --nullomer_filepath --kmer_output_filepath --kmer_length
```
**Main arguments**

###### `--genome_filepath`
Path to the fasta file containing the genome being analyzed. This file must be formatted such that the FASTA headers read ">chr1", ">chr2", "chr3", ...

###### `--nullomer_input_filepath`
Path to the .txt file containing the nullomers absent from the supplied genome.

###### `--mutation_input_filepath`
Path to the file containing the variants being analyzed by this tool.

###### `--mutation_output_filepath`
Path to the output file where the nullomer-causing mutations will be written to.

###### `--nullomer_length`
Length of nullomers being analyzed.

**Optional Arguments**

###### `--chr_col`
Column number of input variant file indicating the variant chromosome (Default=0, for .vcf format)

###### `--pos_col`
Column number of input variant file indicating the variant position (Default=1, for .vcf format)

###### `--ref_col`
Column number of input variant file indicating the reference allele (Default=3, for .vcf format)

###### `alt_col`
Column number of input variant file indicating the alternative allele (Default=4, for .vcf format)

