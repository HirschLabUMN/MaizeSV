import pickle

dirName="/Users/pmonnahan/Documents/Research/MaizeSV/software/gene-key-pipeline/complete_gene_keys/"
outdir="/Users/pmonnahan/Documents/Research/MaizeSV/scripts/accessory/"

def addGene(hom_dict, gene, alt_ref, entry):
	if gene in hom_dict:
		if alt_ref in hom_dict[gene]:
			hom_dict[gene][alt_ref].append(entry)
		else:
			hom_dict[gene][alt_ref] = [entry]
	else:
		hom_dict[gene] = {alt_ref: [entry]}
	return(hom_dict)

def addFile(filename, ref1dict, ref2dict, ref1key, ref2key):
	with open(filename) as r2r:
		for line in r2r:
			line = line.strip().split("\t")
			r1name = line[3]
			r2name = line[7]
			support = line[-1].split(",")
			ref1dict = addGene(ref1dict, r1name, ref2key, [r2name] + support)
			ref2dict = addGene(ref2dict, r2name, ref1key, [r1name] + support)
	return(ref1dict, ref2dict)

def Pickle(Dict, file):
	with open(file, 'wb') as f:
		f.write(pickle.dumps(Dict))

#bhoms = {name: {W22: [gene1,gene2], PHB47: [gene1], PH207: []}}
bhoms = {}
phoms = {}
whoms = {}
Phoms = {}

bhoms, phoms = addFile(dirName + "B73-PH207_gene-key.txt", bhoms, phoms, "B73", "PH207")
bhoms, whoms = addFile(dirName + "B73-W22_GeneKey.txt", bhoms, whoms, "B73", "W22")
bhoms, Phoms = addFile(dirName + "b73-phb47_gene-key.txt", bhoms, Phoms, "B73", "PHB47")
whoms, phoms = addFile(dirName + "w22-ph207_gene-key.txt", whoms, phoms, "W22", "PH207")

# Missing PHB47 to W22 and PHB47 to PH207
Pickle(bhoms, outdir + "b73_homologues_dict.txt")
Pickle(whoms, outdir + "w22_homologues_dict.txt")
Pickle(Phoms, outdir + "phb47_homologues_dict.txt")
Pickle(phoms, outdir + "ph207_homologues_dict.txt")