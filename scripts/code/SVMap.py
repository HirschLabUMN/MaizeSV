

dirName="/Users/pmonnahan/Documents/Research/MaizeSV/software/gene-key-pipeline/complete_gene_keys/"
b2w="B73-W22_GeneKey.txt"
b2P="b73-phb47_gene-key.txt"
w2p="w22-ph207_gene-key.txt"

bgenes={}

bcoords = {chrom: [start1, stop1, name1], [start2, stop2, name2]}
pcoords = {chrom: [start1, stop1, name1]}

bhoms = {name: {W22: [gene1,gene2], PHB47: [gene1], PH207: P[]}}
#bgenes={chrom: [start,stop, name, {W22: [gene1,gene2],PHB47:[gene1],PH207: []}]}
with open(dirName + "B73-PH207_gene-key.txt") as b2p:
	for line in b2p:
		lin = line.split("\t")

		bgenes[line[0]] =

		d
