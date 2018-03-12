import pandas as pd



genes = pd.DataFrame()
with open("/Users/pmonnahan/Documents/Research/Maize/references/Zea_mays.AGPv4.33.genes_only.gff3", 'r') as gff:
	print("start b73")
	for i,line in enumerate(gff):
		if i > 0:
			if i % 1000 == 0: print(i)
			line = line.split()
			genes = genes.append({'ref':'b73', 'name': line[8].split("=")[1].split(";")[0], 'chrom': line[0], 'start': line[3], 'end': line[4]}, ignore_index=True)
print("Number of rows in data frame =", len(genes.index))
with open("/Users/pmonnahan/Documents/Research/Maize/references/ZmaysPH207_443_v1.1.genes_only.gff3", 'r') as gff:
	print("start ph207")
	for i,line in enumerate(gff):
		if i % 1000 == 0: print(i)
		line = line.split()
		genes = genes.append({'ref':'ph207', 'name': line[8].split("=")[-1], 'chrom': line[0], 'start': line[3], 'end': line[4]}, ignore_index=True)
print("Number of rows in data frame =", len(genes.index))
with open("/Users/pmonnahan/Documents/Research/Maize/references/zea_maysw22_core_fixchr.gff", 'r') as gff:
	print("start w22")
	for i,line in enumerate(gff):
		if i % 1000 == 0: print(i)
		line = line.split()
		if line[2] == "gene":
			genes = genes.append({'ref':'w22', 'name': line[8].split("=")[1].split(";")[0], 'chrom': line[0], 'start': line[3], 'end': line[4]}, ignore_index=True)

print("Number of rows in data frame =", len(genes.index))
with open("/Users/pmonnahan/Documents/Research/Maize/references/ZmaysvarPHB47v1.1.gene.gff3", 'r') as gff:
	print("start phb47")
	for i,line in enumerate(gff):
		if i % 1000 == 0: print(i)
		if i > 2:
			line = line.split()
			if line[2] == "gene":
				genes = genes.append({'ref':'phb47', 'name': line[8].split("=")[1].split(";")[0], 'chrom': line[0], 'start': line[3], 'end': line[4]}, ignore_index=True)
print("Number of rows in data frame =", len(genes.index))
genes.to_csv("/Users/pmonnahan/Documents/Research/Maize/AllRefGFF.csv", index = False)