import EnumerateNullomers as en
import itertools

def test_genome_reader():
    # generate an enumerate nullomers object
    enumerator = en.EnumerateNullomers(genome_fasta="sample_data3/tiny_genome.fa", kmer_length=2)
    read_genome = enumerator.genome_sequences
    assert(read_genome == ["ACGTC", "ACTGCA"])

def test_count_2mers():
    # generate an enumerate nullomers object
    enumerator = en.EnumerateNullomers(genome_fasta="sample_data3/tiny_genome.fa", kmer_length=2)
    kmer_counts = enumerator._count_genome_kmers()
    kmers_present = [kmer for kmer in kmer_counts.keys()]
    kmer_counts = [kmer for kmer in kmer_counts.values()]

    # check that the identified present kmers are correct, and also check their counts
    assert(set(kmers_present) == set(["AC", "CG", "GT", "TC", "CT", "TG", "GC", "GA", "AG", "CA"]))
    print(kmer_counts)
    assert(kmer_counts == [3,2,3,1,1,2,2,2,1,1])

def test_2nullomer_id():
    # generate an enumerate nullomers object
    enumerator = en.EnumerateNullomers(genome_fasta="sample_data3/tiny_genome.fa", kmer_length=2)
    nullomers, kmer_counts = enumerator.enumerate()
    assert(nullomers == set(["AT", "AA", "CC", "GG", "TT", "TA"]))

def test_count_3mers():
    # generate an enumerate nullomers object
    enumerator = en.EnumerateNullomers(genome_fasta="sample_data3/tiny_genome.fa", kmer_length=3)
    kmer_counts = enumerator._count_genome_kmers()
    kmers_present = [kmer for kmer in kmer_counts.keys()]
    kmer_counts = [kmer for kmer in kmer_counts.values()]

    assert(set(kmers_present) == set(["ACG", "CGT", "GTC", "ACT", "CTG", "TGC", "GAC", "AGT", "CAG", "GCA"]))
    print(kmer_counts)
    assert(kmer_counts == [2,2,1,1,1,2,2,1,1,1])

def test_3nullomer_id():
    # generate an enumerate nullomers object
    enumerator = en.EnumerateNullomers(genome_fasta="sample_data3/tiny_genome.fa", kmer_length=3)
    nullomers, kmer_counts = enumerator.enumerate()
    kmers_present = [kmer_info[0] for kmer_info in kmer_counts]
    all_3mers = set([''.join(p) for p in itertools.product(["A", "C", "G", "T"], repeat=3)])
    assert(nullomers == all_3mers-set(kmers_present))
    
