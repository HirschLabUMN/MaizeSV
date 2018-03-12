import subprocess
import argparse
import pandas as pd
import re



parser = argparse.ArgumentParser()
parser.add_argument('-vcf', type=str, required=True, help='path to Merged Lumpy VCF, made via lsort and lmerge')
args = parser.parse_args()

jj=0
num_bnd=0
num_other=0
num_del=0
num_dup = 0
num_inv = 0

## Loop Over VCFs ##
with open(args.vcf,'r') as first_vcf:
	for line in first_vcf:
		if line.startswith('##'):
			pass
		elif line.startswith('#'):
			samps = line.split()[9:] # 2 = B73, 23 = ph207, 24 = w22
			# for i,j in enumerate(line.split()):
			# 	print(i,j)
		else:
			line = line.split()
			chrom = line[0]
			pos = int(line[1])
			ID = line[2]
			alt = line[4]
			info = line[7].split(";")
			svtype = info[[i for i, s in enumerate(info) if 'SVTYPE=' in s][0]].split("=")[1] 
			if svtype == "BND":
				kk = re.split('\[|\]',alt)[1].split(":")
				alt_chrom = kk[0]
				alt_pos = kk[1]
			strands = info[1].split("=")[1]
			su = info[[i for i, s in enumerate(info) if 'SU=' in s][0]].split("=")[1]
			pe = info[[i for i, s in enumerate(info) if 'PE=' in s][0]].split("=")[1]
			sr = info[[i for i, s in enumerate(info) if 'SR=' in s][0]].split("=")[1]
			snames = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			svlen = info[[i for i, s in enumerate(info) if 'SNAME=' in s][0]].split("=")[1].split(",")
			num_samps = len(snames)
			if len(ID.split("_")) > 1:
				if ID.split("_")[1]=="1":
					print(num_samps, "BND")
					num_bnd+=1
			else:
				# print(num_samps)
				if svtype == "DEL":
					num_del += 1
					print(num_samps, "DEL")
				elif svtype == "DUP":
					num_dup += 1
					print(num_samps, "DUP")
				elif svtype == "INV":
					num_inv += 1
					print(num_samps, "INV")
				num_other+=1


# print(num_bnd,"BND's", num_bnd/(num_bnd+num_other))
# print(num_del,"DEL's", num_del/(num_bnd+num_other))
# print(num_inv,"INV's", num_inv/(num_bnd+num_other))
# print(num_dup,"DUP's", num_dup/(num_bnd+num_other))
# print(num_other,"Others")

# 			0 SVTYPE=BND
# 1 STRANDS=++:3,--:1
# 2 IMPRECISE
# 3 CIPOS=-19,223
# 4 CIEND=-12,316
# 5 CIPOS95=-10,119
# 6 CIEND95=-4,97
# 7 SU=4
# 8 PE=4
# 9 SR=0
			# if jj < 5:
			# 	for i,j in enumerate(info):
			# 		print(i,j)
			# 		print(su)
			# 	jj+=1
#			get gene name and coords for each reference
			# w22_cmd = 'tabix /Users/pmonnahan/Documents/Research/Maize/W22_Lumpy_Merged.vcf ' + chrom + ":" + str(chrom) + ':' + str(pos)
			
			# if length > 500:
			# 	cnv['dist']=abs(pos - cnv.location)+abs(length - cnv.seq_size)
			# 	min_loc=cnv['dist'].idxmin(axis=1)
			# 	num_deletions += 1
			# 	num_svs +=1
			# 	pp = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)