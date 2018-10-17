
import argparse
import os
from cyvcf2 import VCF
import random
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

def makeDuoPics(vcf_list, sample_list, outdir, bam_dir, samplot_directory, bcftools_executable, num_duos, length_threshold = 100000):
	for i in vcf_list:
		if os.path.exists(i[1]):
			if i[1].endswith("vcf"):
				vcf = VCF(i[1])
				vcf_dir = i[1].split("/")[-1].replace(".vcf","_duos")
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
					else:
						svLen = variant.INFO.get('SVLEN')
					if svLen < length_threshold:
						alts = [j for j, x in enumerate(variant.gt_types) if x == 3]
						refs = [j for j, x in enumerate(variant.gt_types) if x == 0]
						if len(alts) > 2 and len(refs) > 2:
							for k in range(0,num_duos):
								alt = samps[random.choice(alts)]
								ref = samps[random.choice(refs)]
								png_file = f"{svtype}_{variant.CHROM}_{variant.start}_{variant.end}_{alt}_{ref}.png"
								cmd = f"{samplot_directory}/samplot.py -n {alt},{ref} -b {bam_dir}/{alt}{i[0]},{bam_dir}/{ref}{i[0]} -o {Outdir}/{png_file} -s {variant.start} -e {variant.end} -c {variant.CHROM} -a -t {svtype}"
								print(cmd)
			elif i[1].endswith("gz"):
				print("unzip vcf file: ", i[1])		
		else:
			print(i[1], "does not exist") 
	return()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Wrapper for samplot_vcf.sh')
	parser.add_argument('-o', type=str, metavar='output_directory', required=True)
	parser.add_argument('-B', type=str, metavar='bcftools_executable', default="/panfs/roc/msisoft/bcftools/1.6/bin/bcftools")
	parser.add_argument('-S', type=str, metavar='samplot_directory', default="/home/hirschc1/pmonnaha/software/SV-plaudit/Samplot/src/")
	parser.add_argument('-b', type=str, metavar='bam_dir', required=True)
	parser.add_argument('-v', type=str, metavar='vcf_file', required=False, default="-9")
	parser.add_argument('-s', type=str, metavar='bam_suffix', required=False, default="-9", help="sample name in vcf file + suffix  = bamName")
	parser.add_argument('-f', type=str, metavar='vcf_file_list', required=False, default="-9", help="Tab-delimited file with the bam suffix to look for followed by vcf file")
	parser.add_argument('--samps', type=str, metavar='sample_names', required=False, default="-9", help="Sample names can be provided here as a comma-separated list (no spaces).  Otherwise, samples will be retrieved from the header of the vcf")
	parser.add_argument('-d', type=int, metavar='num_duos', required=False, default=-9, help="how many ref/alt sample pairs do you want to generate pictures for")

	args = parser.parse_args()

	if os.path.exists(args.o):
		if args.o.endswith("/") is False: args.o += "/"
	else:
		print(args.o, "does not exist")
	if os.path.exists(args.S):
		if args.S.endswith("/") is False: args.S += "/"
	else:
		print(args.S, "does not exist")

	vcf_list = getVCFlist(args.f, args.v, args.s)
	if args.d != -9:
		makeDuoPics(vcf_list, args.samps, args.o, args.b, args.S, args.B, args.d)



