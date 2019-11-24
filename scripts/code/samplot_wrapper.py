#!/usr/bin/env python
"""This script generates commands for samplot.py, which will generate images of evidence supporting SV calls.  Currently, only variants from Lumpy and Genome STRiP are directly supported, although additional softwares can be easily accommodated.  There is substantial flexibility on the types of images that can be created, but currently it is set up to produce a single picture with the top 3 samples containing the variant and the bottom 3 matching the reference genotype.
	Takes XXX arguments:
	
"""

#Import modules
import argparse
import os
from cyvcf2 import VCF
from numpy import random
import pdb


def getVCFlist(file_list, vcf_file, suffix):
	vcf_list = []
	if (file_list == "-9" and vcf_file == "-9") or (file_list != "-9" and vcf_file != "-9"):
		print("Must provide either a vcf file (-v) and bam suffix (-s) or a file (-f) with paired entries of bam suffix and vcf file")
	elif file_list != "-9":
		assert os.path.exists(file_list), "Could not find file: %r" % file_list
		with open(file_list, 'r') as vcf_file_list:
			for line in vcf_file_list:
				line = line.strip()
				line = line.split()
				assert len(line) == 2, "incorrect format on line: %r" % line
				vcf_list.append(line)
	elif vcf_file != "-9":
		assert suffix != "-9", "Must specify the bam suffix to look for in vcf: %r" % vcf_file
		assert os.path.exists(vcf_file), "Could not find file: %r" % vcf_file
		vcf_list.append([suffix, vcf_file])
	return(vcf_list)

def wholeVCFcommands(vcf_list, sample_list, outdir, bam_dir, samplot_directory, bcftools_executable):
	for i in vcf_list:
		if os.path.exists(i[1]):
			if i[1].endswith("vcf"):
				vcf = VCF(i[1])
				if sample_list == "-9": samps = vcf.samples
				else: samps = sample_list.split(",")
				for sample in samps: 
					vcf_dir = i[1].split('/')[-1].strip('.vcf')
					outdir = f"{outdir}/{vcf_dir}/{sample}"
					bam_file = f"{bam_dir}/{sample}{i[0]}"
					if os.path.exists(args.b + sample + i[0]):
						cmd = f"{samplot_directory}/samplot_vcf.sh -o {outdir} -B {bcftools_executable} -v {i[1]} -S {samplot_directory}/samplot.py {bam_file}"
					else:
						cmd = f"ERROR: Bam file does not exist, {bam_file}"
					try:
						os.makedirs(outdir)
					except FileExistsError:
						cmd = 'WARNING: Output directory already exists;' + cmd
					print(cmd)
			elif i[1].endswith("gz"):
				print("unzip vcf file: ", i[1])		
		else:
			print(i[1], "does not exist") 
	return()

def makeComboPics(vcf_list, sample_list, outdir, bam_dir, samplot_directory, bcftools_executable, num_pics, num_samps, ref_id, length_threshold = 100000):
	for i in vcf_list:
		if os.path.exists(i[1]):
			if i[1].endswith("vcf"):
				suf = i[0]
				vcf = VCF(i[1])
				vcf_dir = i[1].split("/")[-1].replace(".vcf","_combos")
				# pdb.set_trace()
				Outdir = f"{outdir}/{vcf_dir}"
				if not os.path.exists(Outdir): os.mkdir(Outdir)
				if sample_list == "-9": samps = vcf.samples
				else: samps = sample_list.split(",")
				for variant in vcf:
					svtype = variant.INFO.get('SVTYPE')
					if svtype == "CNV":
						svtype = variant.INFO.get('GSCNCATEGORY')
						svLen = variant.INFO.get('GSELENGTH')
						if svtype == "None": print("Change Type to String for GSCNCATEGORY in VCF header")
						genos = variant.format('CN').tolist()
						genos = [x[0] for x in genos]
						if variant.format('FT') is not None:
							filts = [j for j, x in enumerate(variant.format('FT')) if x != "PASS"]
						else: filts = []
						if samps.index(ref_id) in filts: continue
						else: 
							ref_allele = genos[samps.index(ref_id)]
							genos = [0 if x == ref_allele else 3 for x in genos]
							genos = [-9 if j in filts else x for j, x in enumerate(genos)]
					else:
						svLen = variant.INFO.get('SVLEN')
						genos = variant.gt_types
					if svLen < length_threshold:
						alts = [j for j, x in enumerate(genos) if x == 3]
						refs = [j for j, x in enumerate(genos) if x == 0]
						if len(alts) > num_samps and len(refs) > num_samps: #CHANGE NEEDED HERE TO ALLOW FOR 3 AND 3 OR X AND X ALT/REF SAMPS
							for k in range(0, num_pics):
								Samps = [samps[ii] for ii in random.choice(alts, num_samps, replace=False)] + [samps[ii] for ii in random.choice(refs, num_samps, replace=False)]
								# alt = [samps[i] for i in random.choice(alts, num_samps, replace=False)]
								# ref = [samps[i] for i in random.choice(refs, num_samps, replace=False)]
								Bams = [f"{bam_dir}/{x}{suf}" for x in Samps]
								png_file = f"{svtype}_{variant.CHROM}_{variant.start}_{variant.end}.png"
								cmd = f"{samplot_directory}/samplot.py -n {','.join(Samps)} -b {','.join(Bams)} -o {Outdir}/{png_file} -s {variant.start} -e {variant.end} -c {variant.CHROM} -a -t {svtype}"
								print(cmd)
			elif i[1].endswith("gz"):
				print("unzip vcf file: ", i[1])		
		else:
			print(i[1], "does not exist") 
	return()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='This script generates commands for samplot.py, which will generate images of evidence supporting SV calls.  Currently, only variants from Lumpy and Genome STRiP are directly supported, although additional softwares can be easily accommodated.  There is substantial flexibility on the types of images that can be created, but currently it is set up to produce a single picture with the top 3 samples containing the variant and the bottom 3 matching the reference genotype.')
	parser.add_argument('-o', type=str, metavar='output_directory', required=True)
	parser.add_argument('-B', type=str, metavar='bcftools_executable', default="/panfs/roc/msisoft/bcftools/1.6/bin/bcftools")
	parser.add_argument('-S', type=str, metavar='samplot_directory', default="/home/hirschc1/pmonnaha/software/SV-plaudit/Samplot/src/")
	parser.add_argument('-b', type=str, metavar='bam_dir', required=True)
	parser.add_argument('-v', type=str, metavar='vcf_file', required=False, default="-9")
	parser.add_argument('-s', type=str, metavar='bam_suffix', required=False, default="-9", help="sample name in vcf file + suffix  = bamName")
	parser.add_argument('-f', type=str, metavar='vcf_file_list', required=False, default="-9", help="Tab-delimited file with the bam suffix to look for followed by vcf file")
	parser.add_argument('--samps', type=str, metavar='sample_names', required=False, default="-9", help="Sample names can be provided here as a comma-separated list (no spaces).  Otherwise, samples will be retrieved from the header of the vcf")
	parser.add_argument('-N', type=int, metavar='num_pics', required=False, default=1, help="Number of pictures per variant")
	parser.add_argument('-n', type=int, metavar='num_samp', required=False, default=3, help="Number of samples to show for ref and alt alleles")
	parser.add_argument('-c', type=str, metavar='make_combos', required=False, default='true', help="Number of samples to show for ref and alt alleles")
	parser.add_argument('-r', type=str, metavar='reference_id', required=False, default='-9', help="Reference name as it is in VCF.  Used to find reference allele for GSTRiP VCFs.")

	args = parser.parse_args()

	#Check that output directory exists
	if os.path.exists(args.o):
		if args.o.endswith("/") is False: args.o += "/"
	else:
		os.mkdir(args.o)
	#Check that you can find the directory containing the samplot code
	if os.path.exists(args.S):
		if args.S.endswith("/") is False: args.S += "/"
	else:
		print("Unable to find the directory containing samplot code at:", args.S)

	#Ensure proper specification of VCF files
	assert (args.f != "-9" | args.v != "-9"), print("Must provide one (via -v) or more (via -f) VCF files containing SVs that you want to generate pictures from")
	assert (args.f != "-9" & args.v == "-9") | (args.f == "-9" & args.v != "-9"), print("Do not set both -v and -f flags.  Use -v to generate images for a single VCF, and -f to generate images for a list of VCFs specified in a text file.")

	#Parse VCF info to make pictures
	if args.f != "-9":
		vcf_list = getVCFlist(args.f, args.v, args.s)
	elif args.v != "-9": vcf_list = [[args.s, args.v]]

	#Main function to generate commands for generating pictures
	if args.c == 'true':
		makeComboPics(vcf_list, args.samps, args.o, args.b, args.S, args.B, args.N, args.n, args.r)



