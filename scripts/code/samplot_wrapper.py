import vcf
import argparse
import os

parser = argparse.ArgumentParser(description='Wrapper for samplot_vcf.sh')
parser.add_argument('-o', type=str, metavar='output_directory', required=True)
parser.add_argument('-B', type=str, metavar='bcftools_executable', default="/panfs/roc/msisoft/bcftools/1.6/bin/bcftools")
parser.add_argument('-S', type=str, metavar='samplot_directory', default="/home/hirschc1/pmonnaha/software/SV-plaudit/Samplot/src/")
parser.add_argument('-b', type=str, metavar='bam_dir', required=True)
parser.add_argument('-v', type=str, metavar='vcf_file', required=False, default="-9")
parser.add_argument('-s', type=str, metavar='bam_suffix', required=False, default="-9", help="sample name in vcf file + suffix  = bamName")
parser.add_argument('-f', type=str, metavar='vcf_file_list', required=False, default="-9", help="Tab-delimited file with the bam suffix to look for followed by vcf file")

if os.path.exists(args.o):
	if args.o.endswith("/") is False: args.o += "/"
else:
	print(args.o, "does not exist")
if os.path.exists(args.S):
	if args.S.endswith("/") is False: args.S += "/"
else:
	print(args.S, "does not exist")

vcf_list = []
if (args.f == "-9" and args.v == "-9") or (args.f != "-9" and args.v != "-9"):
	print("Must provide either a vcf file (-v) and bam suffix (-s) or a file (-f) with paired entries of bam suffix and vcf file")
elif args.f != "-9":
	if not os.path.exists(args.f):
		print("Could not find file: ", args.f)
		break
	else:
		with open(args.f, 'r') as vcf_file_list:
			for line in vcf_file_list:
				line = line.strip()
				line = line.split("\t")
				assert len(line) == 2
				vcf_list.append(line)
elif args.v != "-9":
	if args.s == "-9":
		print("Must specify the bam suffix to look for vcf: ", args.v)
		break
	if not os.path.exists(args.v):
		print("Could not find file: ", args.v)
		break
	else:
		vcf_list.append([args.v, args.s])

for vcf in vcf_list:
	if os.path.exists(vcf[1]):
		if vcf[1].endswith("vcf"):
			vcf_reader = vcf.Reader(open(vcf[1], 'r'))
			for sample in vcf_reader.samples():
				bams = [s for s in os.listdir(args.b) if sample in s and vcf[0] in s]
				outdir = args.o + vcf[1].strip("vcf") + "/" + sample
				assert len(bams) == 1
				cmd = [args.S + 'samplot_vcf.sh', '-o', outdir, "-B", args.B, "-v", args.v, "-S", args.S + 'samplot.py', bams[0]]
				print(" ".join(cmd))
		elif vcf[1].endswith("gz"):
			print("unzip vcf file: ", vcf[1])		
	else:
		print(vcf[1], "does not exist") 