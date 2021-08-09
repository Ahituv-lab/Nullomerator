import FindAlmostNullomers as fan

def test_generate_left_flank():
    finder = fan.FindAlmostNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", 3)

    assert(finder._generate_left_flank("chr1", 0) == "")
    assert(finder._generate_left_flank("chr1", 1) == "A")
    assert(finder._generate_left_flank("chr1", 2) == "AC")
    assert(finder._generate_left_flank("chr1", 3) == "CG")
    assert(finder._generate_left_flank("chr1", 4) == "GT")

    assert(finder._generate_left_flank("chrX", 0) == "")
    assert(finder._generate_left_flank("chrX", 1) == "A")
    assert(finder._generate_left_flank("chrX", 2) == "AC")
    assert(finder._generate_left_flank("chrX", 3) == "CT")
    assert(finder._generate_left_flank("chrX", 4) == "TG")
    assert(finder._generate_left_flank("chrX", 5) == "GC")


def test_generate_right_flank():
    finder = fan.FindAlmostNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", 3)

    assert(finder._generate_right_flank("chr1", 0, deletion_length=0) == "CG")
    assert(finder._generate_right_flank("chr1", 1, deletion_length=0) == "GT")
    assert(finder._generate_right_flank("chr1", 2, deletion_length=0) == "TC")
    assert(finder._generate_right_flank("chr1", 3, deletion_length=0) == "C")
    assert(finder._generate_right_flank("chr1", 4, deletion_length=0) == "")

    assert(finder._generate_right_flank("chrX", 0, deletion_length=1) == "TG")
    assert(finder._generate_right_flank("chrX", 1, deletion_length=1) == "GC")
    assert(finder._generate_right_flank("chrX", 2, deletion_length=1) == "CA")
    assert(finder._generate_right_flank("chrX", 3, deletion_length=1) == "A")
    assert(finder._generate_right_flank("chrX", 4, deletion_length=1) == "")
    assert(finder._generate_right_flank("chrX", 5, deletion_length=1) == "")


def test_scan_motif_for_nullomers():
    finder = fan.FindAlmostNullomers("test/test_example/tiny_genome.fa", "test/test_example/tiny_genome_3nullomers.txt", 3)

    assert(set(finder._scan_motif_for_nullomers("GACAA")) == set(["ACA", "CAA", "TTG", "TGT"]))
    assert(finder._scan_motif_for_nullomers("ACGTC") == [])

    




#>chr1
#ACGTC

#>chrX
#ACTGCA