import re,os,sys,glob
from Bio import SeqIO
from Bio.Seq import Seq 
import argparse

class ExtractMutationNullomers():
    """
    This class takes in a list of mutations (SNPs, indels, complex subs), and looks to see if any of them generate nullomers given a specified genome and list 
    of nullomers from that genome
    """

    def __init__(self, genome_fasta = "hg38p13_primary.fa", nullomer_file = "All_kmer_words_1_15nt_occurrences_genome_hg38", nullomer_len = 10):
        self.nullomer_set = None
        self.genome_dict = None
        self.genome_length_dict = None
        self.nullomer_len = nullomer_len
        self._genome_reader(genome_fasta)
        self._nullomer_reader(nullomer_file)

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

    def _mutation_reader(self, mutation_file, comments):
        """ 
        reads in a vcf file holding mutation data and returns each line as a list with each word in the line as an element
        """
        datafile=open(mutation_file,"r")
        file_lines = datafile.readlines()
        datafile.close()
        mutation_data = []
        for i in file_lines:
            # if this isnt a comment line, break it up and save it because it contains mutation data
            if i[0] != comments:
                mutation_data+=[i.strip().split('\t')]

        return mutation_data

    def _generate_left_flank(self, chrom_key, pos):
        """
        Captures the genomic sequence preceding a variant.
        """
        # if the mutation position is very close to the beginning or end of the chromosome, the variant motifs need to be handled a bit differently
        # the bases leading up to the mutation site will always be the same, so no conditionals needed regarding changing the left flank sequence
        if pos < self.nullomer_len:
            left_flank = self.genome_dict[chrom_key][:pos - 1]
        else:
            left_flank = self.genome_dict[chrom_key][pos - self.nullomer_len : pos - 1]
        
        return left_flank
    
    def _generate_right_flank(self, chrom_key, pos, ref_len):
        """
        Captures the genomic sequence following a variant.
        """
        # if the ref allele is 1, this is either a snp or an insertion, which are handled the same
        if ref_len == 1:
            # if right flank extends past end of the chromosome, handle this case
            if (pos + self.nullomer_len - 1) >= self.genome_length_dict[chrom_key]:
                right_flank = self.genome_dict[chrom_key][pos:]
            else:
                right_flank = self.genome_dict[chrom_key][pos : pos + self.nullomer_len - 1]
            
        # if ref_len > 1, this is a deletion or complex indel. These cases are handled the same
        elif ref_len > 1:
            if (pos + ref_len + self.nullomer_len - 2) > self.genome_length_dict[chrom_key]:
                right_flank = self.genome_dict[chrom_key][pos + ref_len - 1:]
            else:
                right_flank = self.genome_dict[chrom_key][pos + ref_len - 1: pos + ref_len + self.nullomer_len - 2]
        
        return right_flank

    def scan_mutations(self, mutation_file, chr_col=0, pos_col=1, ref_col=3, alt_col=4, comments="#"):
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
                A list of all of the input variants with relevant information added. For each variant, the chr, pos, reference allele, alternative allele, and nullomers it creates are returned.
            """
        mutation_data = self._mutation_reader(mutation_file, comments="#")
        print("Mutation read-in complete, starting mutation scan")
        nullomer_mutations = []
        for item in mutation_data:
            chrom, pos, ref, alt = item[chr_col], int(item[pos_col]), item[ref_col].upper(), item[alt_col].upper()
            meta_data = [item[i] for i in range(len(item)) if i not in [chr_col, pos_col, ref_col, alt_col]]
            try:
                chrom_int = int(chrom)
                chrom_key = "chr" + chrom
            except: 
                chrom_key = chrom
                
            ref_len = len(ref)
            alt_len = len(alt)
            
            
            # get the left and right flank
            left_flank = self._generate_left_flank(chrom_key, pos)
            right_flank = self._generate_right_flank(chrom_key, pos, ref_len)

            # generate the mutated sequence, and its reverse complement
            if alt == "." or alt == "-":
            	alt = ""
                
            variant_motif = left_flank + alt + right_flank
            variant_rc = str(Seq(variant_motif).reverse_complement())
            motif_len = len(variant_motif)

            # search for nullomers in the new sequence
            found_nullomers = []
            for i in range(0, motif_len-self.nullomer_len+1):
                subsequence = variant_motif[i:i+self.nullomer_len]
                rc_subsequence = variant_rc[i:i+self.nullomer_len]
                if subsequence in self.nullomer_set:
                    found_nullomers.append(subsequence)
                if rc_subsequence in self.nullomer_set:
                    found_nullomers.append(rc_subsequence)
            

            # if at least one new nullomer found, put it with the mutation and append it to a *special* list.
            if pos-self.nullomer_len < 1: # handle edge case where nullomer is right at the beginning of a chromosome
                range_start = 1
                range_end = pos+self.nullomer_len-1
            elif pos+self.nullomer_len > self.genome_length_dict[chrom_key]: # handle edge case where nullomer is right at the end of a chromosome
                range_start = pos-self.nullomer_len+1
                range_end = self.genome_length_dict[chrom_key]
            else:
                range_start = pos-self.nullomer_len
                range_end = pos+self.nullomer_len-1

            if found_nullomers:
                nullomer_mutations.append((chrom, range_start, range_end, ref, alt, variant_motif, found_nullomers))
            else: # if no nullomers are found for that variant, just put "None" in the found_nullomers column
                nullomer_mutations.append((chrom, range_start, range_end, ref, alt, variant_motif, "None"))
            
        return nullomer_mutations
    
    @staticmethod
    def write_output(output_filepath, mutations):
        """
        Writes an output file given a list of lines, where each item in a line is a single item in a list. 
        """
        with open(output_filepath, "w") as mutation_file:
            for mutation in mutations:
                for item in mutation:
                    if isinstance(item, list): 
                        for i,null in enumerate(item):
                            if i == len(item)-1:
                                mutation_file.write(null)
                            else:
                                mutation_file.write(null + ",")
                    elif item == "None":
                        mutation_file.write(item)
                    elif item == "":
                        mutation_file.write("-" + "\t")
                    else:
                        mutation_file.write(str(item) + "\t")
                mutation_file.write("\n")



# Input filename into the object for use
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_file', help='Put in the path to the fasta file containing the genome being analyzed', type=str)
    parser.add_argument('--nullomer_file', help='Put in the path to the .txt file containing the nullomers absent from the supplied genome', type=str)
    parser.add_argument('--mutation_input_file', help='Put in the path to the vcf file containing the mutations being analyzed by this tool', type=str)
    parser.add_argument('--mutation_output_file', help='Put in the path to the output file where the nullomer-causing mutations will be written to', type=str)
    parser.add_argument('--nullomer_length', help="Supply the length of nullomers being analyzed", type=int)

    args = parser.parse_args()
    example_scanner = ExtractMutationNullomers(genome_fasta=args.genome_file, nullomer_file=args.nullomer_file, nullomer_len=args.nullomer_length)

    # change the path to your mutation file. The positional arguments after the path argument mark the columns in which each piece of relevant data lies.
    # for example, looking for mutations stored in a vcf file would look like: 
    #                                           example_scanner.scan_mutations("path_to_vcf", chr_col=0, pos_col=1, ref_col=3, alt_col=4)
    mutations = example_scanner.scan_mutations(args.mutation_input_file, chr_col=0, pos_col=1, ref_col=3, alt_col=4)
    
    # write the output file. if you are using this, just change the output filepath to your desired location 
    ExtractMutationNullomers.write_output(args.mutation_output_file, mutations)
    
            