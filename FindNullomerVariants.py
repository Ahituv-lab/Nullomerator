import re,os,sys,glob,itertools
import numpy as np
import time 
import argparse

class FindNullomerVariants():
    """
    This class finds variants in a provided file that causes nullomers to resurface.
    Takes in a list of population variants that cause nullomers, and a list of mutations that WOULD cause nullomers,
    and finds their intersecting variants
    """

    def __init__(self, population_variants_filepath, nullomer_mutation_filepath, chr_col=0, pos_col=1, ref_col=3, alt_col=4): 
        self.pop_variants = self._mutation_reader(population_variants_filepath, "#", chr_col, pos_col, ref_col, alt_col)
        print("Population variants read-in complete")
        self.nullomer_variants = self._mutation_reader(nullomer_mutation_filepath, nullomers=True)
        print("Nullomer-causing variants read-in complete")


    def _mutation_reader(self, mutation_file, comments="#", chr_col=0, pos_col=1, ref_col=3, alt_col=4,nullomers=False):
        """ 
        reads in a vcf file holding mutation data and returns each line as a list with each word in the line as an element
        """
        datafile=open(mutation_file,"r")
        file_lines = datafile.readlines()
        datafile.close()
        mutation_data = []
        if nullomers==False:
            for i in file_lines:
                # if this isnt a comment line, break it up and save it because it contains mutation data
                if i[0] != comments:
                    mutation_line = i.strip().split('\t')
                    mutation_info = (mutation_line[chr_col].strip("Chrchr"), mutation_line[pos_col], mutation_line[ref_col], mutation_line[alt_col])
                    mutation_data.append(mutation_info)
        else:
            for i in file_lines:
                # if this isnt a comment line, break it up and save it because it contains mutation data
                if i[0] != comments:
                    mutation_line = i.strip().split('\t')
                    mutation_line[0] = mutation_line[0].strip("Chrchr")
                    mutation_data.append(mutation_line)
                    
        return mutation_data
    

    def find_nullomer_variants(self):
        """
        Performs the main operation explained in the class description.
        """
        
        # create a "mutation key" for both of these datasets, of format chr:pos:mutation
        mut_keylist = [line[0] + ":" + line[1] + ":" + line[2] + ":" + line[3]  for line in self.pop_variants]
        nullomer_keylist = [line[0] + ":" + line[1] + ":" + line[3] + ":" + line[4] for line in self.nullomer_variants]
        # use a set for much faster indexing
        mut_variant_keyset = set(mut_keylist)

        true_nullomer_variants = []
        population_nullomer_variants = []
        # go through nullomer keys, and find which ones are in the population mutation bin
        for i, nullomer_key in enumerate(nullomer_keylist):
            if nullomer_key in mut_variant_keyset:
                population_nullomer_variants.append(self.nullomer_variants[i])
            else:
                true_nullomer_variants.append(self.nullomer_variants[i])
        
        return population_nullomer_variants, true_nullomer_variants

def write_output(path, nullomer_variants):
    """
    Function to write a line-by-line output. Handles lists where elements are either strings or sublists
    """
    with open(path, "w") as mutation_file:
        mutation_file.write("#chr\tpos\tid\tref\talt\tvariant_motif\tvariant_motif_interval\tnullomers_created\n")
        for mutation in nullomer_variants:
            for item in mutation:
                if isinstance(item, list): 
                    for i,null in enumerate(item):
                        if i == len(item)-1:
                            mutation_file.write(null)
                        else:
                            mutation_file.write(null + ",")
                elif item == "None":
                    mutation_file.write(item)
                else:
                    mutation_file.write(str(item) + "\t")
            mutation_file.write("\n")


# Input filename into the object for use
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--population_variants_input_filepath', help='Path to the input variant file containing the population variants', type=str)
    parser.add_argument('--nullomer_mutations_input_filepath', help='Path to the input variant file containing mutations that cause neomers', type=str)
    parser.add_argument('--population_variants_output_filepath', help="Path of the output file containing neomer-causing population variants", type=str)
    parser.add_argument('--non_population_variants_output_filepath', help="Path of the output file containing neomer-causing non-population variants", type=str)
    parser.add_argument("--chr_col", help="(Optional) Column number of input variant file indicating the variant chromosome (Default=0, for .vcf format)", default=0, type=int)
    parser.add_argument("--pos_col", help="(Optional) Column number of input variant file indicating the variant position (Default=1, for .vcf format)", default=1,type=int)
    parser.add_argument("--ref_col", help="(Optional) Column number of input variant file indicating the reference allele (Default=3, for .vcf format)", default=3, type=int)
    parser.add_argument("--alt_col", help="(Optional) Column number of input variant file indicating the alternative allele (Default=4, for .vcf format)", default=4, type=int)

    args = parser.parse_args()

    finder = FindNullomerVariants(population_variants_filepath=args.population_variants_input_filepath, 
                                  nullomer_mutation_filepath=args.nullomer_mutations_input_filepath,
                                  chr_col=args.chr_col,
                                  pos_col=args.pos_col,
                                  ref_col=args.ref_col,
                                  alt_col=args.alt_col)
    pop_variants, true_nul_variants = finder.find_nullomer_variants()

    print(len(pop_variants), "of the provided neomer-causing mutations are in the list of population variants")
    write_output(args.population_variants_output_filepath, pop_variants)
    write_output(args.non_population_variants_output_filepath, true_nul_variants)
