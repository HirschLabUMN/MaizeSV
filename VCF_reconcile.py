import subprocess
import argparse
import pandas as pd
import re
import numpy
import gzip



# parser = argparse.ArgumentParser()
# parser.add_argument('-vcf', type=str, required=True, help='path to Merged Lumpy VCF, made via lsort and lmerge')
# args = parser.parse_args()

ph207_vcf = "/Users/pmonnahan/Documents/Research/Maize/PH207_VCFs/PH207_Lumpy_Master.vcf.gz"
b73_vcf = "/Users/pmonnahan/Documents/Research/Maize/B73_VCFs/B73_Lumpy_Master.vcf.gz"
phb47_vcf = "/Users/pmonnahan/Documents/Research/Maize/B73_VCFs/PHB47_PHB47_PHB47.multigt.vcf.gz"
w22_vcf = "/Users/pmonnahan/Documents/Research/Maize/W22_VCFs/W22_Lumpy_Master.vcf.gz"

jj=0
num_bnd=0
num_other=0
num_del=0
num_dup = 0
num_inv = 0
num_A = [0,0,0,0] # These events are where a variant is not within a gene boundary within reference in which it was discovered
num_B = [0,0,0,0] # These events are where a variant is discovered within multiple gene boundaries
ref_discord = [0,0,0,0] # These events are where a variant is heterozygous or alternative homozygous for the reference genotype when mapped to self

out = open("/Users/pmonnahan/Documents/Research/Maize/Lumpy_Consensus.csv", 'w')

# CSV created with make_gff_csv.py
gffs = pd.read_csv("/Users/pmonnahan/Documents/Research/Maize/AllRefGFF.csv")
print(gffs.head)

homs = pd.read_csv("/Users/pmonnahan/Documents/Research/Maize/maize_map.csv")
print(homs.head)

## Loop Over VCFs ##

with gzip.open(b73_vcf,'r') as first_vcf:
	for i,line in enumerate(first_vcf):
		if i % 1000 == 0: print("b73",i)
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

			ref_info = line[11]
			ref_geno = ref_info.split(":")[0]
			if ref_geno != "0/0":
				discord = 1
				ref_discord[0] += 1


			b73_gene1 = gffs.ix[(gffs.ref == 'b73') & (gffs.chrom == chrom) & (gffs.start < pos) & (gffs.end > pos),'name']
			b73_gene2 = pd.DataFrame()
			svtype = info[[i for i, s in enumerate(info) if 'SVTYPE=' in s][0]].split("=")[1] 
			if svtype == "BND":
				kk = re.split('\[|\]',alt)[1].split(":")
				alt_chrom = kk[0]
				alt_pos = int(kk[1])
				b73_gene2 = gffs.ix[(gffs.ref == 'b73') & (gffs.chrom == alt_chrom) & (gffs.start < alt_pos) & (gffs.end > alt_pos),'name']

			strands = info[1].split("=")[1]
			su = info[[i for i, s in enumerate(info) if 'SU=' in s][0]].split("=")[1]
			pe = info[[i for i, s in enumerate(info) if 'PE=' in s][0]].split("=")[1]
			sr = info[[i for i, s in enumerate(info) if 'SR=' in s][0]].split("=")[1]
			snames = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			svlen = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			num_samps = len(snames)

			if len(b73_gene1) == 0 and len(b73_gene2) == 0:
				num_A[0] += 1
			elif len(b73_gene1) != 0 and len(b73_gene2) != 0:
				num_B[0] += 1
			elif len(b73_gene1) == 1 and len(b73_gene2) == 0:
				b73_gene = b73_gene1.iloc[0]
				w22_genes = homs.ix[(homs.gene_id_b73 == b73_gene),'gene_id_w22']
				ph207_genes = homs.ix[(homs.gene_id_b73 == b73_gene),'gene_id_ph207']

				if len(ph207_genes) == 1 and ph207_genes.notnull().iloc[0]:
					ph207_gene = ph207_genes.iloc[0]
					ph207_ginfo = gffs.ix[gffs.name == ph207_gene,]

					tabix_arg = str(ph207_ginfo.chrom.iloc[0]) + ":" + str(ph207_ginfo.start.iloc[0]) + '-' + str(ph207_ginfo.end.iloc[0])
					pp = subprocess.check_output(["tabix", ph207_vcf,tabix_arg])
					# print('pp',pp)

				if len(w22_genes) == 1 and w22_genes.notnull().iloc[0]:
					w22_gene = w22_genes.iloc[0]
					w22_ginfo = gffs.ix[gffs.name == w22_gene,]
					tabix_arg = str(w22_ginfo.chrom.iloc[0]) + ":" + str(w22_ginfo.start.iloc[0]) + '-' + str(w22_ginfo.end.iloc[0])
					ww = subprocess.check_output(["tabix", w22_vcf,tabix_arg])

				if len(ww) > 1:
					num_OL += 1
				if len(pp) > 1:
					num_OL += 1
			# print(num_OL)

with gzip.open(ph207_vcf,'r') as vcf:
	ph207_discord_n = 0
	for i, line in enumerate(vcf):
		if i % 1000 == 0: print("ph207",i)
		num_OL = 0
		ph207_discord = 0
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

			ref_info = line[32]
			ref_geno = ref_info.split(":")[0]
			if ref_geno != "0/0":
				discord = 1
				ref_discord[1] += 1


			ph207_gene1 = gffs.ix[(gffs.ref == 'ph207') & (gffs.chrom == chrom) & (gffs.start < pos) & (gffs.end > pos),'name']
			ph207_gene2 = pd.DataFrame()
			svtype = info[[i for i, s in enumerate(info) if 'SVTYPE=' in s][0]].split("=")[1] 
			if svtype == "BND":
				kk = re.split('\[|\]',alt)[1].split(":")
				alt_chrom = kk[0]
				alt_pos = int(kk[1])
				ph207_gene2 = gffs.ix[(gffs.ref == 'ph207') & (gffs.chrom == alt_chrom) & (gffs.start < alt_pos) & (gffs.end > alt_pos),'name']

			strands = info[1].split("=")[1]
			su = info[[i for i, s in enumerate(info) if 'SU=' in s][0]].split("=")[1]
			pe = info[[i for i, s in enumerate(info) if 'PE=' in s][0]].split("=")[1]
			sr = info[[i for i, s in enumerate(info) if 'SR=' in s][0]].split("=")[1]
			snames = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			svlen = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			num_samps = len(snames)

			if len(ph207_gene1) == 0 and len(ph207_gene2) == 0:
				num_A[1] += 1
			elif len(ph207_gene1) != 0 and len(ph207_gene2) != 0:
				num_B[1] += 1
			elif len(ph207_gene1) == 1 and len(ph207_gene2) == 0:
				ph207_gene = ph207_gene1.iloc[0]
				w22_genes = homs.ix[(homs.gene_id_b73 == ph207_gene),'gene_id_w22']
				b73_genes = homs.ix[(homs.gene_id_b73 == ph207_gene),'gene_id_b73']

				if len(b73_genes) == 1 and b73_genes.notnull().iloc[0]:
					b73_gene = b73_genes.iloc[0]
					b73_ginfo = gffs.ix[gffs.name == b73_gene,]

					tabix_arg = str(b73_ginfo.chrom.iloc[0]) + ":" + str(b73_ginfo.start.iloc[0]) + '-' + str(b73_ginfo.end.iloc[0])

					pp = subprocess.check_output(["tabix", b73_vcf,tabix_arg])
					# print('pp',pp)

				if len(w22_genes) == 1 and w22_genes.notnull().iloc[0]:
					w22_gene = w22_genes.iloc[0]
					w22_ginfo = gffs.ix[gffs.name == w22_gene,]
					tabix_arg = str(w22_ginfo.chrom.iloc[0]) + ":" + str(w22_ginfo.start.iloc[0]) + '-' + str(w22_ginfo.end.iloc[0])
					ww = subprocess.check_output(["tabix", w22_vcf,tabix_arg])

				if len(ww) > 1:
					num_OL += 1
				if len(pp) > 1:
					num_OL += 1

with gzip.open(w22_vcf,'r') as vcf:
	for i,line in vcf:
		if i % 1000 == 0: print("w22",i)
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

			ref_info = line[33]
			ref_geno = ref_info.split(":")[0]
			if ref_geno != "0/0":
				discord = 1
				ref_discord[2] += 1


			b73_gene1 = gffs.ix[(gffs.ref == 'b73') & (gffs.chrom == chrom) & (gffs.start < pos) & (gffs.end > pos),'name']
			b73_gene2 = pd.DataFrame()
			svtype = info[[i for i, s in enumerate(info) if 'SVTYPE=' in s][0]].split("=")[1] 
			if svtype == "BND":
				kk = re.split('\[|\]',alt)[1].split(":")
				alt_chrom = kk[0]
				alt_pos = int(kk[1])
				b73_gene2 = gffs.ix[(gffs.ref == 'b73') & (gffs.chrom == alt_chrom) & (gffs.start < alt_pos) & (gffs.end > alt_pos),'name']

			strands = info[1].split("=")[1]
			su = info[[i for i, s in enumerate(info) if 'SU=' in s][0]].split("=")[1]
			pe = info[[i for i, s in enumerate(info) if 'PE=' in s][0]].split("=")[1]
			sr = info[[i for i, s in enumerate(info) if 'SR=' in s][0]].split("=")[1]
			snames = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			svlen = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			num_samps = len(snames)

			if len(b73_gene1) == 0 and len(b73_gene2) == 0:
				num_A += 1
			elif len(b73_gene1) != 0 and len(b73_gene2) != 0:
				num_B += 1
			elif len(b73_gene1) == 1 and len(b73_gene2) == 0:
				b73_gene = b73_gene1.iloc[0]
				w22_genes = homs.ix[(homs.gene_id_b73 == b73_gene),'gene_id_w22']
				ph207_genes = homs.ix[(homs.gene_id_b73 == b73_gene),'gene_id_ph207']

				if len(ph207_genes) == 1 and ph207_genes.notnull().iloc[0]:
					ph207_gene = ph207_genes.iloc[0]
					ph207_ginfo = gffs.ix[gffs.name == ph207_gene,]

					if chrom != "10":
						tabix_arg = 'chr0' + chrom + ":" + str(ph207_ginfo.start.iloc[0]) + '-' + str(ph207_ginfo.end.iloc[0])
					else:
						tabix_arg = 'chr' + chrom + ":" + str(ph207_ginfo.start.iloc[0]) + '-' + str(ph207_ginfo.end.iloc[0])
					pp = subprocess.check_output(["tabix", ph207_vcf,tabix_arg])
					# print('pp',pp)

				if len(w22_genes) == 1 and w22_genes.notnull().iloc[0]:
					w22_gene = w22_genes.iloc[0]
					w22_ginfo = gffs.ix[gffs.name == w22_gene,]
					tabix_arg = 'chr' + chrom + ":" + str(w22_ginfo.start.iloc[0]) + '-' + str(w22_ginfo.end.iloc[0])
					ww = subprocess.check_output(["tabix", w22_vcf,tabix_arg])

				if len(ww) > 1:
					num_OL += 1
				if len(pp) > 1:
					num_OL += 1


			
print("b73 discordant sites =", b73_discord_n)