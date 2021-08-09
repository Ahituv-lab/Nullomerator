import FindNullomerVariants as fnv

def test_mutation_reader():
    finder = fnv.FindNullomerVariants("test/test_example/tiny_genome_pop_variants.txt", "test/test_example/nullomer-generating_mutations.tsv")
    known_pop_variants = [("1","1","A","C"),
                            ("1","1","A","AG"),
                            ("1","1","AC","A"),
                            ("1","3","G","T"),
                            ("1","4","TC","T"),
                            ("1","5","C","A"),
                            ("1","5","C","T"),
                            ("1","5","C","CA"),
                            ("X","1","A","G"),
                            ("X","3","T","A"),
                            ("X","6","A","T")]
    assert(finder.pop_variants == known_pop_variants)

    ex_nul_variant = finder.nullomer_variants[0]
    assert(ex_nul_variant==["1", "1", ".", "A", "G", "GCG", "GCG,CGC"])
    assert(len(finder.nullomer_variants) == 82)

def test_nullomer_finder():
    finder = fnv.FindNullomerVariants("test/test_example/tiny_genome_pop_variants.txt", "test/test_example/nullomer-generating_mutations.tsv")
    pop_variants, true_nul_variants = finder.find_nullomer_variants()

    assert(len(true_nul_variants) == 73)
    assert(len(pop_variants) == 9)


#chr	pos	id	ref	alt
#1	1	.	A	C
#1	1	.	A	AG
#1	1	.	AC	A
#1	3	.	G	T
#1	4	.	TC	T
#1	5	.	C	A
#1	5	.	C	T
#1	5	.	C	CA
#X	1	.	A	G
#X	3	.	T	A
#X	6	.	A	T

