python3 ../EnumerateNullomers.py\
 --genome_filepath tiny_genome.fa\
 --nullomer_output_filepath tiny_genome_4nullomers.txt\
 --kmer_output_filepath tiny_genome_4mers.txt\
 --kmer_length 4

python3 ../ExtractMutationNullomers.py\
 --genome_filepath tiny_genome.fa\
 --nullomer_input_filepath tiny_genome_4nullomers.txt\
 --mutation_input_filepath tiny_genome_patient_variants.txt\
 --mutation_output_filepath tiny_genome_patient_neomer_variants.txt\
 --nullomer_length 4

python3 ../FindAlmostNullomers.py\
 --genome_filepath tiny_genome.fa \
 --nullomer_input_filepath tiny_genome_4nullomers.txt\
 --output_filepath tiny_genome_all_neomer_variants.txt\
 --nullomer_length 4

python3 ../FindNullomerVariants.py\
 --population_variants_input_filepath tiny_genome_pop_variants.txt\
 --nullomer_mutations_input_filepath tiny_genome_all_neomer_variants.txt\
 --population_variants_output_filepath tiny_genome_neomer_pop_variants.txt\
 --non_population_variants_output_filepath tiny_genome_neomer_variants.txt