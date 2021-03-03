import re,os,sys,glob
from Bio import SeqIO
from Bio.Seq import Seq 

class ExtractMutationNullomers():

    def __init__(self, genome_fasta = "hg38p13_primary.fa", nullomer_file = "All_kmer_words_1_15nt_occurrences_genome_hg38", nullomer_len = 0):
        self.nullomer_set = None
        self.genome_dict = None
        self.max_null_length = None
        self.nullomer_len = nullomer_len
        self._genome_reader(genome_fasta)
        self._nullomer_reader(nullomer_file)

    def _genome_reader(self, genome_fasta):
        """
        reads in a genome in fasta format, and sets a dictionary where each key is a chrom and each entry is a string of the sequence

        Arguments:
            genome_fasta::str
                Path to a genome fasta file that only holds the 24 primary chromosomes, with fasta headings ">chr1", ">chr2",...">chrX", "chrY"
        
        Returns:
            None
        """
        genome_dict = {}
        for section in SeqIO.parse(genome_fasta, "fasta"):
            genome_dict[str(section.description)] = str(section.seq)
        self.genome_dict = genome_dict
        print("Genome read-in complete")
        

    def _nullomer_reader(self, nullomer_file):
        """
        Reads in a two-column file, with kmers in one col and the number of times they occur in the genome in the second col. Sets a set of kmers with count == 0.

        Arguments:
            nullomer_file::str
                Path to a file of kmers, where the left column is the kmer and the right column is the count of its occurrences in the genome
        
        Returns:
            None
        """
        nullomer_set = set()
        max_null = 0
        with open(nullomer_file) as infile:
            # look for nullomers of any length
            if self.nullomer_len == 0:
                for line in infile:
                    kmer, count = line.split()
                    kmer_len = len(kmer)
                    if int(count)==0:
                    # if the occurence of a kmer is 0, it is a nullomer and add it to the set
                        nullomer_set.add(kmer)
                        if kmer_len > max_null: 
                            self.max_null_length = kmer_len
            else: # look for nullomers of only a certain length
                self.max_null_length = self.nullomer_len
                for line in infile:    
                    kmer, count = line.split()
                    kmer_len = len(kmer)
                    if int(count) == 0 and kmer_len == self.nullomer_len:
                        nullomer_set.add(kmer)

            self.nullomer_set = nullomer_set
            print("Nullomer read-in complete")

    def _mutation_reader(self, mutation_file):
        """ 
        reads in a vcf file holding mutation data and returns each line as a list with each word in the line as an element
        """
        datafile=open(mutation_file,"r")
        file_lines = datafile.readlines()
        datafile.close()
        mutation_data = []
        for i in file_lines[1:]:
            items = [i.strip().split('\t')]
            mutation_data+=[i.strip().split('\t')]     
        return mutation_data


    def scan_mutations(self, mutation_file, chr_col=0, pos_col=1, ref_col=3, alt_col=4):
        """
        Takes in a vcf file of mutations, and looks for whether the create a nullomer in the genome

        Arguments:
            mutation_file::str
                A path to a vcf file storing mutation data
            chr_col::int
                Column of the vcf where chr data is stored. Default=0
            pos_col::int
                Column of the vcf where position data is stored. Default=1
            ref_col::int
                Column of the vcf where the reference allele is stored. Default=3
            alt_col::int
                Column of the vcf where the alternative allele is stored. Default=4
        
        Returns:
            nullomer_mutations::list
                A list of variants that create nullomers. For each variant, the chr, pos, reference allele, alternative allele, and nullomers it creates are returned.
            """
        mutation_data = self._mutation_reader(mutation_file)
        print("Mutation read-in complete, starting mutation scan")
        nullomer_mutations = []
        for item in mutation_data:
            chrom, pos, ref, alt = item[chr_col], int(item[pos_col]), item[ref_col], item[alt_col]
            cancer_type, patient_id = item[0], item[1]
            chrom_key = "chr" + chrom
            ref_len = len(ref)
            alt_len = len(alt)
            
            # the bases leading up to the mutation site will always be the same, so no conditionals needed
            left_flank = self.genome_dict[chrom_key][pos - self.max_null_length : pos - 1].upper()

            # if the ref allele is 1, this is either a snp or an insertion, which are handled the same
            if ref_len == 1:
                right_flank = self.genome_dict[chrom_key][pos : pos + self.max_null_length - 1].upper()
            
            # if ref_len > 1, this is a deletion or complex indel. These cases are handled the same
            elif ref_len > 1:
                right_flank = self.genome_dict[chrom_key][pos + ref_len - 1: pos + ref_len + self.max_null_length - 2].upper()

            # generate the mutated sequence, and its reverse complement
            variant_motif = left_flank + alt + right_flank
            variant_rc = str(Seq(variant_motif).reverse_complement())
            motif_len = len(variant_motif)

            # search for nullomers in the new sequence
            found_nullomers = []
            if self.nullomer_len == 0:
                for seq_len in range(1, self.max_null_length + 1):
                    for i in range(0, motif_len-seq_len+1):
                        subsequence = variant_motif[i:i+seq_len]
                        rc_subsequence = variant_rc[i:i+seq_len]

                        if subsequence in self.nullomer_set:
                            found_nullomers.append(subsequence)
                        if rc_subsequence in self.nullomer_set:
                            found_nullomers.append(rc_subsequence)
                        
            # if the nullomer length is fixed, the computation is slightly different
            else:
                for i in range(0, motif_len-self.nullomer_len):
                    subsequence = variant_motif[i:i+self.nullomer_len]
                    rc_subsequence = variant_rc[i:i+self.nullomer_len]
                    if subsequence in self.nullomer_set:
                        found_nullomers.append(subsequence)
                    if rc_subsequence in self.nullomer_set:
                        found_nullomers.append(rc_subsequence)

            # if at least one new nullomer found, put it with the mutation and append it to a *special* list
            if len(found_nullomers) != 0:
                # get length of longest nullomer
                max_len = 0
                for null in found_nullomers:
                    length = len(null)
                    if  length > max_len:
                        max_len = length
                
                # each line should be of format: chrom, start, end, cancer type patient id REF ALT seq nullomers
                nullomer_mutations.append((chrom, pos-max_len, pos+max_len, cancer_type, patient_id, ref, alt, variant_motif, found_nullomers))

        return nullomer_mutations



# Input filename into the object for use
if __name__ == "__main__":
    # change the genome_fasta path to your genome file
    example_scanner = ExtractMutationNullomers(genome_fasta="genome_files/chr1_hg37.fa", nullomer_file="All_kmer_words_1_15nt_occurrences_genome_hg38", nullomer_len=0)
    # change the path to your mutation file
    mutations = example_scanner.scan_mutations("genome_files/mutation_example.txt", chr_col=2, pos_col=3, ref_col=4, alt_col=5)
    
    # write the output file
    with open("mutation_nullomers.txt", "w") as mutation_file:
        for mutation in mutations:
            for item in mutation:
                if isinstance(item, list):
                    for null in item:
                        mutation_file.write(null + ",")
                else:
                    mutation_file.write(str(item) + "\t")
            mutation_file.write("\n")
            