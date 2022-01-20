.. find_almost_nullomers:
  
====================
FindAlmostNullomers
====================

|

**Identifies all possible single base pair substitutions, insertions, or deletions, genome-wide that cause the resurfacing of at least one nullomer.**

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
   
   * - --nullomer_output_filepath
     - Path to the output .txt file where the nullomers absent from the supplied genome will be written, with one nullomer written on each line.

   * - ---nullomer_length
     - Length of nullomers being analyzed.

|

|

**Outputs:**

* VCF-compatible file containing all variants genome-wide that cause at least one nullomer to resurface. 

|
