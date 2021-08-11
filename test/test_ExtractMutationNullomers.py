import ExtractMutationNullomers as emn
import itertools

def test_genome_reader():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_2nullomers.txt", nullomer_len=2)
    chromosome_seqs = [val for val in extractor.genome_dict.values()]
    chromosome_labels = [val for val in extractor.genome_dict.keys()]
    chromosome_lengths = [val for val in extractor.genome_length_dict.values()]

    assert(chromosome_seqs == ["ACGTC", "ACTGCA"])
    assert(chromosome_labels == ["chr1", "chrX"])
    assert(chromosome_lengths == [5,6])

def test_nullomer_reader():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_2nullomers.txt", nullomer_len=2)
    assert(extractor.nullomer_set == set(["AT", "AA", "CC", "GG", "TT", "TA"]))

def test_mutation_reader():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_2nullomers.txt", nullomer_len=2)
    mutations = extractor._mutation_reader("test/test_example/tiny_genome_patient_variants.txt", comments="#")

    assert(len(mutations)==14)

    mut_chrs = [mut[0] for mut in mutations]
    assert(mut_chrs == ["1","1","1","1","1","1","1","1","1","1","1","X","X","X"])

def test_generate_left_flank():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", nullomer_len=3)
    mutations = extractor._mutation_reader("test/test_example/tiny_genome_patient_variants.txt", comments="#")

    generated_left_flanks = []
    for item in mutations:
        chrom, pos, ref, alt = item[0], int(item[1]), item[3].upper(), item[4].upper()
        chrom_key = "chr" + chrom
        ref_len = len(ref)
        alt_len = len(alt)

        lf = extractor._generate_left_flank(chrom_key,pos, ref_len,alt_len)
        generated_left_flanks.append(lf)
    assert(generated_left_flanks == ["", "", "","", "AC","AC", "C", "C", "G", "GT", "T", "","C","GC"])

def test_generate_right_flank():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", nullomer_len=3)
    mutations = extractor._mutation_reader("test/test_example/tiny_genome_patient_variants.txt", comments="#")

    generated_right_flanks = []
    for item in mutations:
        chrom, pos, ref, alt = item[0], int(item[1]), item[3].upper(), item[4].upper()
        chrom_key = "chr" + chrom
        ref_len = len(ref)

        lf = extractor._generate_right_flank(chrom_key,pos, ref_len)
        generated_right_flanks.append(lf)
    
    assert(generated_right_flanks == ["CG", "CG", "GT","GT", "C","TC", "TC", "C", "", "", "", "CT","A", ""])



def test_nullomer_identification():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", nullomer_len=3)
    mutations = extractor.scan_mutations("test/test_example/tiny_genome_patient_variants.txt", chr_col=0, pos_col=1, ref_col=3, alt_col=4)

    variant_motifs=[]
    nullomers_created = []
    for mut in mutations:
        variant_motifs.append(mut[-3])
        nullomers_created.append(mut[-1])

    assert(variant_motifs == ["CCG", "AGCG", "AGT", "TGGT","ACGACC", "ACTTC", "CGTTC", "CGC", "GT", "GTA", "TCA", "CCT", "CTA", "GCT"])

def test_variant_motif_range():
    extractor = emn.ExtractMutationNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", nullomer_len=3)
    mutations = extractor.scan_mutations("test/test_example/tiny_genome_patient_variants.txt", chr_col=0, pos_col=1, ref_col=2, alt_col=3)

    var_ranges = [mut[-2] for mut in mutations]
    
    assert(var_ranges == ["1:3", "1:3", "1:3","1:3","1:5", "1:5", "1:5", "1:5", "2:5", "3:5", "3:5", "1:3", "1:5", "4:6"])

    
