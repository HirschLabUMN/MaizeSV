import pickle
import argparse
import os
import pdb


def addGene(hom_dict, gene, alt_ref, entry):
	if gene == "Zm00014a000149": pdb.set_trace()
	if gene in hom_dict:
		if alt_ref in hom_dict[gene] and entry not in hom_dict[gene][alt_ref]:
			hom_dict[gene][alt_ref].append(entry)
		else:
			hom_dict[gene][alt_ref] = [entry]
	else:
		hom_dict[gene] = {alt_ref: [entry]}
	return(hom_dict)

def addFile(in_file, ref1dict, ref2dict, ref1key, ref2key):

    with open(in_file, 'r') as inFile:
        for i, line in enumerate(inFile):
            if i > 0 and line[0] != "#":
                line = line.split()
                genes = line[11].split(";")
                for gene in genes:
                	gene = gene.split(",")
                	gene = [x.replace('"','') for x in gene]
                	ref2dict = addGene(ref2dict, gene[0], ref1key, line[0])
                	ref1dict = addGene(ref1dict, line[0], ref2key, gene[0])
    return(ref1dict, ref2dict)

def Pickle(Dict, file):
	with open(file, 'wb') as f:
		f.write(pickle.dumps(Dict))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='geneKey_fileList', required=True, help='tab delimited file with ref names and gene key file')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True)
    # parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')
    args = parser.parse_args()

    if not os.path.exists(args.o): os.mkdir(args.o)

    HOMS = {}
    FILES = {}
    with open(args.f, 'r') as key_file:
        for line in key_file:
            ref1, ref2, file = line.strip().split()
            FILES[file] = [ref1, ref2]

    
    for file in FILES:
    	ref1 = FILES[file][0]
    	ref2 = FILES[file][1]
    	try:
    		HOMS[ref1]
    	except KeyError:
    		HOMS[ref1] = {}
    	try:
    		HOMS[ref2]
    	except KeyError:
    		HOMS[ref2] = {}
    	HOMS[ref1], HOMS[ref2] = addFile(file, HOMS[ref1], HOMS[ref2], ref1, ref2)
    	# HOMS[ref2] = addFile(file, HOMS[ref2], HOMS[ref1], ref2, ref1)

    for ref in HOMS:
    	Pickle(HOMS[ref], f"{args.o}/{ref}_homologues_dict.txt")


#bhoms = {name: {W22: [gene1,gene2], PHB47: [gene1], PH207: []}}
# bhoms = {}
# phoms = {}
# whoms = {}
# Phoms = {}

# bhoms, phoms = addFile(dirName + "B73-PH207_gene-key.txt", bhoms, phoms, "B73", "PH207")
# bhoms, whoms = addFile(dirName + "B73-W22_GeneKey.txt", bhoms, whoms, "B73", "W22")
# bhoms, Phoms = addFile(dirName + "b73-phb47_gene-key.txt", bhoms, Phoms, "B73", "PHB47")
# whoms, phoms = addFile(dirName + "w22-ph207_gene-key.txt", whoms, phoms, "W22", "PH207")

# # Missing PHB47 to W22 and PHB47 to PH207
# Pickle(bhoms, outdir + "b73_homologues_dict.txt")
# Pickle(whoms, outdir + "w22_homologues_dict.txt")
# Pickle(Phoms, outdir + "phb47_homologues_dict.txt")
# Pickle(phoms, outdir + "ph207_homologues_dict.txt")

