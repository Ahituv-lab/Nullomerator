.. find_nullomer_variants:
  
====================
FindNullomerVariants
====================

|

**Removes nullomers that are likely to result from common variants in a user specified variant VCF file.**

|

|

.. list-table:: **Required arguments:**
   :header-rows: 1

   * - Argument
     - Explanation
     

   * - --population_variants_input_filepath
     - Path to the input variant file containing common population variants.

   * - --nullomer_mutations_input_filepath
     - Path to the input variant file containing mutations that cause neomers.
   
   * - --population_variants_output_filepath
     - Path of the output file containing neomer-causing population variants.
   
   * - --non_population_variants_output_filepath
     - Path of the output file containing nullomer-causing non-population variants.

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

   * - --freq_col   
     - Column number of input population variant file indicating frequency of the alternative allele (MAF) (Default=5)

|

|

**Outputs:**

* VCF-compatible file containing all variants from `nullomer_mutations_input_filepath` that do not contain a common popoulation variant. 

|