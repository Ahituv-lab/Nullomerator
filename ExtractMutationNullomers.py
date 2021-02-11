import time
import os
from Bio import SeqIO
"""
This program uses a reference genome and a list of nullomers to obtain a list of all possible mutations of the genome 
that could cause the nullomer to resurface. 
NOTE: Certain nullomers might resurface in the same location from two different mutations. For example, for the genome
AGCA, the nullomer AGA can resurface from the substitution of C by A, the deletion of C and the insertion of an A 
before C. This program will display all three mutations as possible causes of this nullomer.
"""
N = 12  # DEFINE NULLOMER LENGTH HERE
nul_file_path = "genome_not_nullomers/12bp_nullomers/All_kmer_words_1_15nt_occurrences_genome_hg38_12bps"
genome_file_path = "chr11.fa"
# If output_tab_separated == false, the output will be formatted as an aligned table for easier human reading.
# If true it will be a tab separated file for easier machine reading
output_tab_separated = True

def mutate(motif, position):
    """This function returns all possible mutations of its input
    Args:
        motif(string): motif of length n+1, where n is the nullomer
        position(int): specifies the position of the first base in motif
    Returns:
        mutations (set): A set of tuples. Each tuple contains the position of the mutation in the chromosome/genome,
        the mutation, the resulting kmer and the position of the mutation in the resulting kmer. Position might be
        needed to be adapted to a tuple to also include chromosome.
        The first position of the genome and of a kmer is 0. The position of an insertion indicates the position of the
        base after the insertion.
        NEW - the tuple  contains the position of the mutated base within the mutation.
    """
    mutations = set()
    n = len(motif)
    for i in range(1, n - 1):
        # Add all possible substitutions
        for base in {"A", "C", "G", "T"} - {motif[i]}:
            mutations.add(
                (position + i - 1, motif[i - 1:i + 2] + ":" + base, motif[1:i] + base + motif[i + 1:-1], i - 1))
        # Add all possible insertions
        for base in {"A", "C", "G", "T"}:
            mutations.add(
                (position + i - 1, motif[i - 1] + "-" + motif[i] + ":" + base, motif[1:i] + base + motif[i:-2], i - 1))
        # Add all possible deletions
        mutations.add((position + i - 1, motif[i - 1:i + 2] + ":-", motif[1:i] + motif[i + 1:], i - 1))
        # BUG TO FIX: At the last position of the genome (or chromosome) the deletion will give the wrong sequence
        # Adding an extra if to handle just this case seems computationally wasteful - maybe need to split before
        # calling mutate and handle last motif separately
    return mutations

def nul_file_reader(path):
    """
    Read the file in the path given, assuming the first entry in each column is a nullomer and returns
    a set with all nullomers
    """
    with open(path) as f1:
        nullomers = set()
        for line in f1:
            nullomers.add(line.split()[0])
    return nullomers

def analysis(name, sequence):
    seq_len = len(sequence)
    mut_nul_set = set()
    for i in range(1, len(sequence) - N):
        mut_set = mutate(sequence[i - 1:i + N + 1], i - 1)
        for mutation in mut_set:
            if mutation[2] in nul_set:
                mutation = (name,) + mutation
                mut_nul_set.add(mutation)
        if i % 100000 == 0:
            print("Currently testing for mutations in position {} out of {} in chromosome {}".format(i, seq_len, name))
    return mut_nul_set


# NEEDS TO BE ADAPTED to the fasta format
def fasta_reader(path):

    return fasta_sequences


def output_writer(mut_nul_set, time_stamp):
    """
    This function takes a set of mutations and writes it in tab separated format in the output file.

    :param: mutations: A set containing tuples, each tuple representing a mutation causing a nullomer to resurfrace.
    :param: time_stamp: Time stamp to append to the output file - needs to be the same for all iterations of the
    function
    The function is agnostic to the format of the tuple. Currently the tuple has the form as described in mutate()
    :return: None
    """

    file_name = "output2/" + time_stamp + "mut_to_nul" + str(N) + "_out" ".txt"
    if not os.path.isfile(file_name):
        with open(file_name, "x") as f2:
            # Add header to output file
            header = ("Chromosome", "Position in genome", "Mutation", "Found nullomer", "Position of mutation within "
                                                                                        "nullomer")
            # I know we'll change the format back, I just wanted my output to look organised
            if output_tab_separated:
                f2.write("\t".join(header)+"\n")
            else:
                f2.write('{:<20} {:<20} {:<20} {:<20} {:<20}'.format(*header) + "\n")

    with open(file_name, "a") as f2:
        for mutation in mut_nul_set:
            if output_tab_separated:
                f2.write('\t'.join([str(x) for x in mutation])+"\n")
            else:
                f2.write('{:<20} {:<20} {:<20} {:<20} {:<20}'.format(*mutation) + "\n")


if __name__ == "__main__":
    time_stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())  # To add timestamp to output file names
    if not os.path.exists("output2/"):
        os.makedirs("output2/")

    # Load the nullomers in a python set
    tic = time.perf_counter()
    nul_set = nul_file_reader(nul_file_path)
    toc = time.perf_counter()

    # Measure time it took to load nullomers
    nul_read_time = "Loaded the nullomer file to memory in {:5f} seconds\n".format(toc - tic)

    file_name = "output2/" + time_stamp + "_mut_to_nul" + str(N) + "_timer" + ".txt"
    with open(file_name, "a") as time_file:
        time_file.write(nul_read_time)

    # Read the genome file
    fasta_sequences = SeqIO.parse(open(genome_file_path), 'fasta')

    # Run the main analysis function and store the output in a file
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        # Add a dummy character at the start and end of the sequence, since mutate discards the two exterior
        # characters
        sequence = "s" + sequence + "s"

        tic = time.perf_counter()
        mut_nul_set = analysis(name, sequence)
        toc = time.perf_counter()

        output_writer(mut_nul_set, time_stamp)

        # Measure the time it took to run the analysis function on the entire sequence
        mut_analysis_time = "Iterated through all possible mutations of sequence \"{}\" in {:0.5f} seconds\n".format(name, toc - tic)

        # Write the time it took to process the sequence in a file
        file_name = "output2/" + time_stamp + "_mut_to_nul" + str(N) + "_timer" + ".txt"
        with open(file_name, "a") as time_file:
            time_file.write(mut_analysis_time)
