import subprocess
import argparse
import pandas as pd
import re
import numpy
import gzip



# parser = argparse.ArgumentParser()
# parser.add_argument('-vcf', type=str, required=True, help='path to Merged Lumpy VCF, made via lsort and lmerge')
# args = parser.parse_args()


vcfs = {'b73':["/Users/pmonnahan/Documents/Research/Maize/B73_VCFs/B73_Lumpy_Master.vcf.gz",2],'ph207':["/Users/pmonnahan/Documents/Research/Maize/PH207_VCFs/PH207_Lumpy_Master.vcf.gz",23],'w22':["/Users/pmonnahan/Documents/Research/Maize/W22_VCFs/W22_Lumpy_Master.vcf.gz",24]} # The number indicates the genotype index for each reference

jj=0
num_bnd=0
num_other=0
num_del=0
num_dup = 0
num_inv = 0
num_A = {'b73':0,'ph207':0,'w22':0,'phb47':0} # These events are where a variant is not within a gene boundary within reference in which it was discovered
num_B = {'b73':0,'ph207':0,'w22':0,'phb47':0} # These events are where a variant is discovered within multiple gene boundaries
ref_discord = {'b73':0,'ph207':0,'w22':0,'phb47':0} # These events are where a variant is heterozygous or alternative homozygous for the reference genotype when mapped to self

out = open("/Users/pmonnahan/Documents/Research/Maize/Lumpy_Consensus.csv", 'w')

# CSV created with make_gff_csv.py
gffs = pd.read_csv("/Users/pmonnahan/Documents/Research/Maize/AllRefGFF.csv")
# gffs = gffs.set_index('name').T.to_dict('list')
# print(gffs['Zm00008a029620'])
homs = pd.read_csv("/Users/pmonnahan/Documents/Research/Maize/maize_map.csv")

## Loop Over VCFs ##
for name, values in vcfs.items():
	alt_vcfs = {x: vcfs[x] for x in vcfs if x not in [name,'phb47']}
	t_gffs = gffs.ix[(gffs.ref == name)]
	with gzip.open(values[0],'r') as vcf:
		old_chrom = -99
		for i,line in enumerate(vcf):
			if i % 1000 == 0: print(name,i)
			num_OL = 0
			discord = 0
			line = line.decode('utf-8')
			if line.startswith('##'):
				pass
			elif line.startswith('#'):
				samps = line.split()[9:] # 2 = B73, 23 = ph207, 24 = w22
				# for i,j in enumerate(line.split()):
				# 	print(i,j)
			else:
				# print(line)
				line = line.split()
				chrom = line[0]
				pos = int(line[1])
				ID = line[2]
				alt = line[4]
				info = line[7].split(";")


				ref_info = line[values[1] + 11] # 11 is offset needed to start at samples genotypes in vcf
				ref_geno = ref_info.split(":")[0]
				if ref_geno != "0/0":
					discord = 1
					ref_discord[name] += 1
				else:
					if chrom != old_chrom:
						gene1 = t_gffs.ix[(t_gffs.chrom == chrom)]
						
					svtype = info[[i for i, s in enumerate(info) if 'SVTYPE=' in s][0]].split("=")[1] 

					if svtype == "BND":
						kk = re.split('\[|\]',alt)[1].split(":")
						alt_chrom = kk[0]
						alt_pos = int(kk[1])
					else:
						svlen = info[[i for i, s in enumerate(info) if 'SVLEN=' in s][0]].split("=")[1]
						alt_chrom = chrom
						alt_pos = pos + abs(int(svlen)) 
					
					gene1 = t_gffs.ix[(t_gffs.chrom == chrom) & (t_gffs.start < pos) & (t_gffs.end > pos),'name']
					gene2 = t_gffs.ix[(t_gffs.chrom == alt_chrom) & (t_gffs.start < alt_pos) & (t_gffs.end > alt_pos),'name']
					if gene1.equals(gene2) is True:
						gene2 = []

					strands = info[1].split("=")[1]
					su = info[[i for i, s in enumerate(info) if 'SU=' in s][0]].split("=")[1]
					pe = info[[i for i, s in enumerate(info) if 'PE=' in s][0]].split("=")[1]
					sr = info[[i for i, s in enumerate(info) if 'SR=' in s][0]].split("=")[1]
					snames = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
					num_samps = len(snames)

					if len(gene1) == 0 and len(gene2) == 0:
						num_A[name] += 1
					elif len(gene1) != 0 and len(gene2) != 0:
						num_B[name] += 1
					elif len(gene1) == 1 and len(gene2) == 0:
						gene = gene1.iloc[0]
						alt1_genes = homs.ix[(homs['gene_id_' + name] == gene),'gene_id_' + list(alt_vcfs.keys())[0]]
						alt2_genes = homs.ix[(homs['gene_id_' + name] == gene),'gene_id_' + list(alt_vcfs.keys())[1]]
						pp = []
						if len(alt1_genes) == 1 and alt1_genes.notnull().iloc[0]:
							alt1_gene = alt1_genes.iloc[0]
							alt1_ginfo = gffs.ix[gffs.name == alt1_gene,]

							tabix_arg = str(alt1_ginfo.chrom.iloc[0]) + ":" + str(alt1_ginfo.start.iloc[0]) + '-' + str(alt1_ginfo.end.iloc[0])
							pp = subprocess.check_output(["tabix", alt_vcfs[list(alt_vcfs.keys())[0]][0],tabix_arg])

							# print('pp',pp)
						ww = []
						if len(alt2_genes) == 1 and alt2_genes.notnull().iloc[0]:
							alt2_gene = alt2_genes.iloc[0]
							alt2_ginfo = gffs.ix[gffs.name == alt2_gene,]
							tabix_arg = str(alt2_ginfo.chrom.iloc[0]) + ":" + str(alt2_ginfo.start.iloc[0]) + '-' + str(alt2_ginfo.end.iloc[0])
							ww = subprocess.check_output(["tabix", alt_vcfs[list(alt_vcfs.keys())[1]][0],tabix_arg])

						if len(ww) > 1:
							num_OL += 1
						if len(pp) > 1:
							num_OL += 1
				# print(num_OL)


			
print("b73 discordant sites =", ref_discord['b73'])
print("ph207 discordant sites =", ref_discord['ph2077'])
print("phb47 discordant sites =", ref_discord['phb47'])
print("w22 discordant sites =", ref_discord['w22'])