import re,os,sys,glob,itertools
import numpy as np
import time 

class FindNullomerVariants():
    """
    This class finds variants in a provided file that causes nullomers to resurface.
    Takes in a list of population variants that cause nullomers, and a list of mutations that WOULD cause nullomers,
    and finds their intersecting variants
    """

    def __init__(self): pass

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

    def find_nullomer_variants(self, mutation_filepath, nullomer_mutation_filepath):
        """
        Performs the main operation explained in the class description.
        Takes in a filepath to population variants in "mutation_filepath", and a list of mutations that are 
        known to cause nullomers("nullomer_mutation_filepath")
        """

        muts = self._mutation_reader(mutation_filepath, comments="#")
        nullomer_muts = self._mutation_reader(nullomer_mutation_filepath, comments="#")

        # make X and Y chromosomes "23" and "24" for pop variants
        for i in range(len(muts)): 
            # change X to "23" and Y to "24"
            if muts[i][0] == "X":
                nullomer_muts[i,0] = "23"
            elif muts[i][0] == "Y":
                nullomer_muts[i][0] = "24"
        
        # remove "chr" from the chromosome labels for each observation
        for i in range(len(nullomer_muts)): 
            nullomer_muts[i][0] = nullomer_muts[i][0].strip("chrChr")
            # change X to "23" and Y to "24"
            if nullomer_muts[i][0] == "X":
                nullomer_muts[i][0] = "23"
            elif nullomer_muts[i][0] == "Y":
                nullomer_muts[i][0] = "24"
        


        # create a "mutation key" for both of these datasets, of format chr:pos:mutation
        mut_keylist = [line[0] + ":" + line[1] + ":" + line[2] for line in muts]
        nullomer_keylist = [line[0] + ":" + line[1] + ":" + line[2] for line in nullomer_muts]
        # use a set for way faster indexing
        nullomer_keyset = set(nullomer_keylist)

        nullomer_variants = []
        # go through mutation keys, and find which ones are in the nullomer-causing mutation bin
        # this is a really slow way to do this, but since we need information from both datasets, manual index tracking is necessary
        for i, mut_key in enumerate(mut_keylist):
            if mut_key in nullomer_keyset:
                # if we find a hit, find it's corresponding nullomer data
                nullomer_infos = []
                for j, nul_key in enumerate(nullomer_keylist):
                    if mut_key == nul_key: 
                        nullomer_infos.append((nullomer_muts[j][3], nullomer_muts[j][4]))

                variant_line = (muts[i][0], muts[i][1], muts[i][2], muts[i][4], [tup[0] for tup in nullomer_infos], [tup[1] for tup in nullomer_infos])
                nullomer_variants.append(variant_line)
        
        return nullomer_variants

def write_output(path, nullomer_variants):
    """
    Function to write a line-by-line output. Handles lists where elements are either strings or sublists
    """
    with open(path, "w") as mutation_file:
    # write header line
        mutation_file.write("#chr"+"\t"+"pos"+"\t"+"mutation"+"\t"+"nullomer_index"+"\t"+"nullomer(s)"+"\t"+"variant_pos_in_nullomer\n")
        for line in nullomer_variants:
            line_string = ""
            for i, item in enumerate(line):
                # if a list, write each item of the list in csv format in one column
                if isinstance(item, list):
                    list_str = ""
                    for sub_item in item:
                        list_str+=(sub_item+",")
                    list_str = list_str.strip(",")
                    line_string += list_str+"\t"
                # if a regular item, just write the item
                else:
                    line_string += item+"\t"
            line_string = line_string.strip()
            mutation_file.write(line_string+"\n")


# Input filename into the object for use
if __name__ == "__main__":
    example = FindNullomerVariants()
    variants = example.find_nullomer_variants(mutation_filepath="nullomers13mutpos.tsv.100k_lines", nullomer_mutation_filepath="nullomer_causing_mutations_head.tsv")
    print(str(len(variants)) + " nullomer-causing variants re-identified")

    write_output("example_output_Py.txt", variants)
