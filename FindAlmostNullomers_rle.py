import re,os,sys,glob,itertools
from Bio.Seq import Seq 
from Bio import SeqIO
import time 
import argparse

class FindAlmostNullomers():
    """
    This class uses a reference genome and a list of nullomers to obtain a list of all possible mutations of the genome 
    that could cause the nullomer to resurface. 
    """

    def __init__(self, genome_file, nullomer_file, nullomer_length):
        self.nullomer_len = nullomer_length
        self._nullomer_reader(nullomer_file)
        self._genome_reader(genome_file)

    def _genome_reader(self, genome_fasta):
        """
        reads in a genome in fasta format, and sets a dictionary where each key is a chrom and each entry is a string of the sequence.
        Also makes a dictionary that holds the lengths of each chromosome

        Arguments:
            genome_fasta::str
                Path to a genome fasta file that only holds the 24 primary chromosomes, with fasta headings ">chr1", ">chr2",...">chrX", "chrY"
        
        Returns:
            None
        """
        genome_dict = {}
        genome_length_dict = {}
        for section in SeqIO.parse(genome_fasta, "fasta"):
            section_description = str(section.description).split()
            section_header = section_description[0]
            chrom_sequence = str(section.seq).upper()
            genome_dict[section_header] = chrom_sequence
            genome_length_dict[section_header] = len(chrom_sequence)
        self.genome_dict = genome_dict
        self.genome_length_dict = genome_length_dict
        print("Genome read-in complete")
    
    def _nullomer_reader(self, nullomer_file):
        """
        Reads in a one-column file, with each line being a nullomer string

        Arguments:
            nullomer_file::str
                Path to a file of nullomer, where each line is a nullomer
        
        Returns:
            None
        """
        nullomer_set = set()
        with open(nullomer_file) as infile:
            for nullomer in infile:
                nullomer = nullomer.strip()
                nullomer_set.add(nullomer)
        self.nullomer_set = nullomer_set
        print("Nullomer read-in complete")
    
    @staticmethod
    def check_bases(seq):
        """
        Checks whether all of the bases in seq are one of the four canonical bases. Returns True if this is case, else False. 
        """
        for base in seq:
            if base != "A" and base != "C" and base != "G" and base != "T" and base != "":
                return False
        return True

    
    def _generate_left_flank(self, chrom_key, pos):
        """
        Captures the genomic sequence preceding a variant.
        """
        # if the mutation position is very close to the beginning or end of the chromosome, the variant motifs need to be handled a bit differently
        # the bases leading up to the mutation site will always be the same, so no conditionals needed regarding changing the left flank sequence
        if pos+1 < self.nullomer_len:
            left_flank = self.genome_dict[chrom_key][:pos]
        else:
            left_flank = self.genome_dict[chrom_key][pos - self.nullomer_len + 1: pos]
        
        return left_flank
    
    def _generate_right_flank(self, chrom_key, pos, deletion_length=0):
        """
        Captures the genomic sequence following a variant.
        """
        # if ref_len > 1, this is a deletion or complex indel. These cases are handled the same
        if (pos + deletion_length + self.nullomer_len) > self.genome_length_dict[chrom_key]:
            right_flank = self.genome_dict[chrom_key][pos+deletion_length+1:]
        else:
            right_flank = self.genome_dict[chrom_key][pos+deletion_length+1: pos+deletion_length+self.nullomer_len]

        return right_flank

    def _scan_motif_for_nullomers(self, motif):
        """
        given a list of DNA motifs, this function finds the ones that has nullomers in them. 

        Arguments:
            motif::string
                Variant motif that is being scanned for nullomers
            
        Returns:
            nullomers_created::[string]
                List of nullomers that are contained in the variant motif and the variant motif's reverse complement.
        """
        nullomers_created = []
        motif_len = len(motif)
        # make sure to check reverse complements as well
        motif_rc = str(Seq(motif).reverse_complement())

        # scan the motifs for nullomers
        for i in range(0, motif_len-self.nullomer_len+1):
            subsequence = motif[i:i+self.nullomer_len]
            rc_subsequence = motif_rc[i:i+self.nullomer_len]
            # save nullomer-causing mutations in "chr, pos, ref, alt, nullomers_created" format
            if subsequence in self.nullomer_set:
                nullomers_created.append(subsequence)
            if rc_subsequence in self.nullomer_set:
                nullomers_created.append(rc_subsequence)
        
        return nullomers_created

    
    def _generate_snp_insertion_variant_sequences_and_scan(self, chrom_number, pos, left_flank, ref_base, right_flank, all_nullomer_variants_list):
        """
        """
        
        dna_bases = {"A", "C", "T", "G"}
        # generate variant motifs that occur from SNPs and small (1 bp) insertions. Then scan these motifs for nullomers.
        for dna_base in dna_bases:
            #generate the variant motif and scan it for nullomers
            # handle snps
            if dna_base != ref_base:
                snp_variant_motif = left_flank + dna_base + right_flank
                snp_nullomers = self._scan_motif_for_nullomers(snp_variant_motif)
                if snp_nullomers:
                    # put together list of needed information (chrom, pos, ref, alt, variant_motif, nullomers_created)
                    info = (chrom_number, pos+1, ref_base, dna_base, snp_variant_motif, snp_nullomers)
                    all_nullomer_variants_list.append(info)
            
            # handle insertions
            insertion_variant_motif = left_flank + ref_base + dna_base + right_flank
            insertion_nullomers = self._scan_motif_for_nullomers(insertion_variant_motif)
            if insertion_nullomers:
                # put together list of needed information (chrom, pos, ref, alt, variant_motif, nullomers_created)
                info = (chrom_number, pos+1, ref_base, ref_base+dna_base, insertion_variant_motif, insertion_nullomers)
                all_nullomer_variants_list.append(info)
    
    def _generate_deletion_variant_sequences_and_scan(self,chrom_number, pos, left_flank, ref_genotype, right_flank, all_nullomer_variants_list):
        # generate variant motif from small (1 bp) deletion. Add this to nullomer causing variant list if nullomers are created
        ref_base = ref_genotype[0]
        deletion_variant_motif = left_flank + ref_base + right_flank
        deletion_nullomers = self._scan_motif_for_nullomers(deletion_variant_motif)
        if deletion_nullomers:
            # put together list of needed information (chrom, pos, ref, alt, variant_motif, nullomers_created)
            # handle the alt allele differently if this is the last position on the chromosome
            info = (chrom_number, pos+1, ref_genotype, ref_base, deletion_variant_motif, deletion_nullomers)
            all_nullomer_variants_list.append(info)



    def find_almost_nullomers(self, indel_size=1):
        """
        With a nullomer set and genome, this method generates a list of mutations that cause a nullomer to arise in the given genome
        """
        all_nullomer_variants = []
        dna_bases = {"A", "C", "G", "T"}
        
        for chromosome, sequence in self.genome_dict.items():
            chrom_number = chromosome.strip("chr")
            
            # loop through each base in the chromosome/sequence 
            for pos, base in enumerate(sequence):
                
                left_flank = self._generate_left_flank(chromosome, pos)
                # if this flank contains a non-canonical base, skip onto the next base
                if not self.check_bases(left_flank):
                    continue

                # handle snps and insertions
                snp_insertion_right_flank = self._generate_right_flank(chromosome, pos)
                if self.check_bases(snp_insertion_right_flank):
                    self._generate_snp_insertion_variant_sequences_and_scan(chrom_number, pos, left_flank, base, snp_insertion_right_flank,all_nullomer_variants)
                
                # handle deletions
                deletion_right_flank = self._generate_right_flank(chromosome, pos, deletion_length=indel_size)
                if deletion_right_flank and self.check_bases(deletion_right_flank):
                    deletion_ref_genotype = sequence[pos:pos+indel_size+1]
                    self._generate_deletion_variant_sequences_and_scan(chrom_number, pos, left_flank, deletion_ref_genotype, deletion_right_flank, all_nullomer_variants)
            
        return all_nullomer_variants

    @staticmethod
    def write_output(output_filepath, mutations):
        """
        Writes an output file given a list of lines, where each item in a line is a single item in a list. 
        """
        with open(output_filepath, "w") as mutation_file:
            mutation_file.write("#chr\tpos\tid\tref\talt\tvariant_motif\tnullomers_created\n")
            for mutation in mutations:
                chr, pos, ref, alt, variant_motif, nullomers_created = mutation
                mutation_file.write(str(chr)+"\t"+str(pos)+"\t"+".\t"+ref+"\t"+alt+"\t"+variant_motif+"\t")
                for i,null in enumerate(nullomers_created):
                    if i == len(nullomers_created)-1:
                        mutation_file.write(null)
                    else:
                        mutation_file.write(null + ",")
                mutation_file.write("\n")

            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_file', help='Put in the path to the fasta file containing the genome being analyzed', type=str)
    parser.add_argument('--nullomer_file', help='Put in the path to the .txt file containing the nullomers absent from the supplied genome', type=str)
    parser.add_argument('--nullomer_length', help="Supply the length of nullomer being analyzed", type=int)

    args = parser.parse_args()

    nullomerator = FindAlmostNullomers(genome_file=args.genome_file, 
                                       nullomer_file=args.nullomer_file,
                                       nullomer_length=args.nullomer_length)
    found_nullomers = nullomerator.find_almost_nullomers()
    FindAlmostNullomers.write_output("sample_data3/nullomer-generating_mutations.tsv", found_nullomers)