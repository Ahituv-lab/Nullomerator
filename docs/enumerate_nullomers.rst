.. enumerate_nullomers:
  
====================
EnumerateNullomers
====================

|

**Extracts all nullomers, and enumerates all present kmers, of a specified kmer length from a fasta sample.**


|

|

.. list-table:: **Required arguments:**
   :header-rows: 1

   * - Argument
     - Explanation
     

   * - --genome_filepath
     - Path to the fasta file containing the genome being analyzed.\ :sup:`*` 

   * - --nullomer_output_filepath
     - Path to the output .txt file where the nullomers absent from the supplied genome will be written, with one nullomer written on each line.
   
   * - --kmer_output_filepath
     - Path to the output .tsv file where the kmers, and their corresponding occurrence counts, from the supplied genome will be written.
   
   * - --kmer_length
     - Length of kmers/nullomers to be enumerated.

|

|

**Outputs:**

* Text file containing a list of nullomers absent from the supplied genome file.
* Text file containing every kmer occuring at least once in the genome file supplied, accompanied by a count of how often each occur

|

.. note::

  \ :sup:`*` \ This file must be formatted such that the FASTA headers read ">chr1", ">chr2", ">chr3", ...

| 