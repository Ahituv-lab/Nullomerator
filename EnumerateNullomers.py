import re,os,sys,glob,itertools
from Bio.Seq import Seq 
from Bio import SeqIO
from collections import Counter
import time 

class EnumerateNullomers():

    def __init__(self, genome_fasta, kmer_length, no_trivial = False, range = True):
        self.kmer_len = kmer_length
        self.no_trivial = no_trivial
        self.range = range
        self._genome_reader(genome_fasta)
    
    @staticmethod
    def check_bases(seq):
        """
        Checks whether all of the bases in seq are one of the four canonical bases. Returns True if this is case, else False. 
        """
        for base in seq:
            if base != "A" and base != "C" and base != "G" and base != "T":
                return False
        return True

    
    def _genome_reader(self, genome_fasta):
        """
        reads in a genome in fasta format, and sets a dictionary where each key is a chrom and each entry is a string of the sequence

        Arguments:
             genome_fasta::str
                Path to a genome fasta file that only holds the 24 primary chromosomes, with fasta headings ">chr1", ">chr2",...">chrX", "chrY"
            
        Returns:
            None
        """
        self.genome_sequences = []
        for section in SeqIO.parse(genome_fasta, "fasta"):
            # strip sequences of unknown bases and make everything uppercase
            self.genome_sequences.append(str(section.seq).upper())
        print("Genome read-in complete")

    def _count_genome_kmers(self):
        """
        Generates a list of all the kmers contained in the provided genome. each kmer is accompanied by a count of how many times it appears.
        Also considers the other strand of the genome by generating counts from the reverse complement as well

        Arguments:
            None
        
        Returns:
            None
        """
        # generate reverse complements
        genome_rcs = [str(Seq(sequence).reverse_complement()) for sequence in self.genome_sequences]
        kmer_counts = Counter()

        # count kmer occurrences for a range of kmer lengths from 1...self.kmer_len
        if self.range == True:
            for kmer_len in range(1, self.kmer_len + 1):
                for seq in self.genome_sequences:
                    for i in range(0, len(seq)+1-kmer_len):
                        motif = seq[i:i+kmer_len]
                        if self.check_bases(motif):
                            kmer_counts[motif] += 1
                for seq in genome_rcs:
                    for im in range(0, len(seq)+1-kmer_len):
                        motif = seq[im:im+kmer_len]
                        if self.check_bases(motif):
                            kmer_counts[motif] += 1
        # count kmer occurrences for only kmers of length self.kmer_len
        else:
            for seq in self.genome_sequences:
                for i in range(0, len(seq)+1-self.kmer_len):
                    motif = seq[i:i+self.kmer_len]
                    if self.check_bases(motif):
                        kmer_counts[motif] += 1
            for seq in genome_rcs:
                for im in range(0, len(seq)+1-self.kmer_len):
                    motif = seq[im:im+self.kmer_len]
                    if self.check_bases(motif):
                        kmer_counts[motif] += 1

        return kmer_counts
    
    def _generate_all_kmers(self):
        """
        Generates every single possible kmer of length self.kmer_len

        Arguments:
            None
        
        Returns:
            all_kmers::[str]

        """
        bases = ["A", "C", "T", "G"]
        all_kmers = []
        # only generate kmers of certain length if range == False
        if self.range == True:
            for kmer_len in range(1, self.kmer_len+1):
                all_kmers+=[''.join(p) for p in itertools.product(bases, repeat=kmer_len)]
        else:
            all_kmers+=[''.join(p) for p in itertools.product(bases, repeat=self.kmer_len)]

        return all_kmers

    def enumerate(self):
        """
        Wrapper function for the main functionality of the class, which is to count the instances of every kmer that appears in the 
        provided genome file, and to determine which kmers of a specified length never appear in the genome file

        Arguments:
            None

        Returns:
            nullomers:[str]
                List of nullomers for the given genome file
            genome_kmers:[array]
                List of tuples, where each tuple is a given kmer with the count of how many times it appears in the genome file
        """
        # get the kmer counts for the genome
        genome_kmers = self._count_genome_kmers()
        # generate every possible kmer of a certain length
        possible_kmers = self._generate_all_kmers()

        # find which kmers do not appear in the genome
        nullomers = set(possible_kmers) - set(genome_kmers.elements()) 
        genome_kmers = [(kmer, genome_kmers[kmer]) for kmer in genome_kmers]

        return nullomers, genome_kmers


if __name__ == "__main__":
    start = time.time()
    # if youre using this script, you can put your desired kmer length below
    nully = EnumerateNullomers(genome_fasta = "genome_files/test.fa", kmer_length = 3)
    nullomers, kmers = nully.enumerate()

    # write output file of genome kmers (kmer with count on same line)
    with open("nullomer+kmer_files/test_kmers.txt", "w") as kmer_file:
        for kmer, count in kmers:
            kmer_file.write(kmer + "\t" + str(count) + "\n")

    # write output file of nullomers (one nullomer per line)
    with open("nullomer+kmer_files/test_nullomers.txt", "w") as nullomer_file:
        for nullomer in nullomers:
            nullomer_file.write(nullomer + "\n")

    end = time.time()
    print((end-start)/60, "minutes elapsed")