#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script prepares the gene-keys produced by JM's gene-key pipeline (https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline), so that they can be used to link SVs across references (SVMap.py)
 Takes the following arguments:
    -f :  REQUIRED: tab delimited file 3 columns: reference1 reference2 ref1-ref1_gene-key-file.   Note that there is a directionality to the gene-key files.  E.g. The gene-key file, 500kb_W22_B73_AllbyAll_res.txt, specifically corresponds to W22 as reference1 and B73 as reference2
    -o :  REQUIRED: Full path to output directory in which the pickled VCFs will be written

"""

#Import necessary modules
import pickle
import argparse
import os
import pdb


def addGene(hom_dict, gene, alt_ref, entry):
    '''This helper function adds homolog data to the correct location as specified in the arguments'''
    if gene in hom_dict: #check if an entry for this gene already exists
        if alt_ref in hom_dict[gene] and entry not in hom_dict[gene][alt_ref]: #Check if there is an entry for the alt_ref and NO entry for current gene
            hom_dict[gene][alt_ref].append(entry) #If so, add current entry to list of genes from this alt_ref that is homologous to gene.  Thus, multiple genes from alt match this gene in ref.
        else:
            hom_dict[gene][alt_ref] = [entry]
    else: #Create new key:item pair if no entry yet exists for this gene
        hom_dict[gene] = {alt_ref: [entry]}
    return(hom_dict)

def addFile(in_file, ref1dict, ref2dict, ref1key, ref2key):
    '''For a particular gene-key file, this function parses the data and adds it to appropriate entries in the HOMS dictionary'''
    with open(in_file, 'r') as inFile:
        for i, line in enumerate(inFile):
            if i > 0 and line[0] != "#": #Skip header lines in gene-key fiile
                line = line.split()
                ref_gene = line[0]
                alt_genes = line[11].split(";") #Grab all of the genes in the alt annotation corresponding this ref_gene entry.
                for alt_gene in alt_genes: #Loop over homologs and add to correct location in HOMS
                    alt_gene = alt_gene.split(",")
                    alt_gene = [x.replace('"','') for x in alt_gene]
                    if alt_gene[0] == '0' or ref_gene == '0': #No matches 
                        continue
                    else:
                       ref2dict = addGene(ref2dict, alt_gene, ref1key, ref_gene) #returns a modified ref2dict which now includes data passed in arguments
                       ref1dict = addGene(ref1dict, ref_gene, ref2key, alt_gene)
    return(ref1dict, ref2dict)

def Pickle(Dict, file):
	with open(file, 'wb') as f:
		f.write(pickle.dumps(Dict))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script prepares the gene-keys produced by JM's gene-key pipeline (https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline), so that they can be used to link SVs across references (SVMap.py)")
    parser.add_argument('-f', type=str, metavar='geneKey_fileList', required=True, help='tab delimited file 3 columns: reference1 reference2 ref1-ref1_gene-key-file.   Note that there is a directionality to the gene-key files.  E.g. The gene-key file, 500kb_W22_B73_AllbyAll_res.txt, specifically corresponds to W22 as reference1 and B73 as reference2')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.o): os.mkdir(args.o) #Make output directory if it does not exist

    #Large dictionary of dictionaries (of dictionaries) with first index corresponding to a reference genotype, secondary indexes corresponding to gene names from the reference annotation, and tertiary indexes corresponding to the alternative references.
    #HOMS[ref][gene][alt_ref] will be a list of gene IDs from the alt_ref annotation that correspond the the 'gene' in the current annotation.
    HOMS = {} 

    #Dictionary holding the reference genomes that correspond to each file.
    FILES = {}
    with open(args.f, 'r') as key_file:
        for line in key_file:
            ref1, ref2, file = line.strip().split()
            FILES[file] = [ref1, ref2]

    
    #Loop over all gene-key files and add entries to appropriate location in HOMS dictionary.
    for file in FILES:
    	ref1 = FILES[file][0] #Get references associated with this file
    	ref2 = FILES[file][1]
    	try: #See if entry exists in HOMS for this reference, and if not, create one.
    		HOMS[ref1]
    	except KeyError:
    		HOMS[ref1] = {}
    	try: #See if entry exists in HOMS for the alt reference, and if not, create one.
    		HOMS[ref2]
    	except KeyError:
    		HOMS[ref2] = {}
    	HOMS[ref1], HOMS[ref2] = addFile(file, HOMS[ref1], HOMS[ref2], ref1, ref2) #Main function that adds info from the file to the HOMS entries for each reference.

    #After data from all files has been stored in HOMS, write the homolog data for each reference to a compressed (pickled) file for fast retrieval within SVMap.py
    for ref in HOMS:
    	Pickle(HOMS[ref], f"{args.o}/{ref}_homologues_dict.txt")

