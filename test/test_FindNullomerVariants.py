from re import L
import FindNullomerVariants as fnv

def test_mutation_reader():
    finder = fnv.FindNullomerVariants("test/test_example/tiny_genome_pop_variants.txt", "test/test_example/nullomer-generating_mutations.tsv")
    known_pop_variants = [("1","1","A","C", ".3"),
                            ("1","1","A","AG", ".25"),
                            ("1","1","AC","A", ".1"),
                            ("1","3","G","T", ".2"),
                            ("1","4","TC","T", ".02"),
                            ("1","5","C","A", ".3"),
                            ("1","5","C","T", ".45"),
                            ("1","5","C","CA", ".12"),
                            ("X","1","A","G", ".07"),
                            ("X","3","T","A", ".05"),
                            ("X","6","A","T", ".25")]
    assert(finder.pop_variants == known_pop_variants)

    ex_nul_variant = finder.nullomer_variants[0]
    assert(ex_nul_variant==["1", "1", ".", "A", "G", "GCG", "GCG,CGC"])
    assert(len(finder.nullomer_variants) == 82)

def test_nullomer_finder():
    finder = fnv.FindNullomerVariants("test/test_example/tiny_genome_pop_variants.txt", "test/test_example/nullomer-generating_mutations.tsv")
    pop_variants, true_nul_variants = finder.find_nullomer_variants()

    assert(len(true_nul_variants) == 73)
    assert(len(pop_variants) == 9)

def test_frequency_calculator():
    finder = fnv.FindNullomerVariants("test/test_example/tiny_genome_pop_variants.txt", "test/test_example/nullomer-generating_mutations.tsv")
    pop_variants, true_nul_variants = finder.find_nullomer_variants()
    pop_nullomer_freq_dict = finder.calculate_neomer_probabilities(pop_variants)
    known_freqs = {
        "AAG":.2,
        "AAC":.45,
        "AGC":.50303,
        "ACA":.05,
        "TTC":.2,
        "TAC":.3,
        "TGA":.12,
        "TCA":.12,
        "TGT":.05,
        "CGC":.25,
        "CCG":.3,
        "CGG":.3,
        "CTT":.2,
        "GCG":.25,
        "GCT":.50303,
        "GAA":.2,
        "GTT":.45,
        "GTA":.3

    }
    for null in pop_nullomer_freq_dict.keys():
        known = known_freqs[null]
        test = round(1-pop_nullomer_freq_dict[null],5)
        assert(known == test)


