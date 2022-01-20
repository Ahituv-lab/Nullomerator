.. extract_mutation_nullomers:
  
========================
ExtractMutationNullomers
========================

|

**Finds all mutations that cause the resurfacing of a list of nullomers.**

|

|

.. list-table:: **Required arguments:**
   :header-rows: 1

   * - Argument
     - Explanation
     

   * - --genome_filepath
     - Path to the fasta file containing the genome being analyzed.\ :sup:`*` 

   * - --nullomer_input_filepath
     - Path to the .txt file containing the nullomers absent from the supplied genome.
   
   * - --mutation_input_filepath
     - Path to the file containing variants to be analyzed.
   
   * - --mutation_output_filepath
     - Path to the output file where the neomer-causing mutations will be written to.

   * - --nullomer_length
     - Length of nullomers being analyzed.

|

| 

.. list-table:: **Optional arguments:**
   :header-rows: 1

   * - Argument
     - Explanation

   * - --chr_col
     - Column number of input variant file indicating the variant chromosome (Default=0, for .vcf format)

   * - --pos_col
     - Column number of input variant file indicating the variant position (Default=1, for .vcf format)

   * - --ref_col
     - Column number of input variant file indicating the reference allele (Default=3, for .vcf format)

   * - --alt_col   
     - Column number of input variant file indicating the alternative allele (Default=4, for .vcf format)

|

|

**Outputs:**

* VCF-compatible file containing variants from provided variant file that cause at least one nullomer to resurface. 

|


.. note::

  \ :sup:`*` \ This file must be formatted such that the FASTA headers read ">chr1", ">chr2", ">chr3", ...

| 