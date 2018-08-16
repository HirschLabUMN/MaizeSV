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
parser.add_argument('--samps', type=str, metavar='sample_names', required=False, default="-9", help="Sample names can be provided here as a comma-separated list (no spaces).  Otherwise, samples will be retrieved from the header of the vcf")

args = parser.parse_args()

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
	assert os.path.exists(args.f), "Could not find file: %r" % args.f
	with open(args.f, 'r') as vcf_file_list:
		for line in vcf_file_list:
			line = line.strip()
			line = line.split()
			assert len(line) == 2, "incorrect format on line: %r" % line
			vcf_list.append(line)
elif args.v != "-9":
	assert args.s != "-9", "Must specify the bam suffix to look for in vcf: %r" % args.v
	assert os.path.exists(args.v), "Could not find file: %r" % args.v
	vcf_list.append([args.s, args.v])


for i in vcf_list:
	if os.path.exists(i[1]):
		if i[1].endswith("vcf"):
			if args.samps == "-9":
				vcf_reader = vcf.Reader(open(i[1], 'r'))
				samps = vcf_reader.samples
			else:
				samps = args.samps.split(",")
			for sample in samps: 
				outdir = args.o + i[1].split("/")[-1].strip(".vcf") + "/" + sample
				if os.path.exists(args.b + sample + i[0]):
					cmd = [args.S + 'samplot_vcf.sh', '-o', outdir, "-B", args.B, "-v", i[1], "-S", args.S + 'samplot.py', args.b + sample + i[0]]
				else:
					cmd = ["ERROR: Bam file does not exist,", args.b + sample + i[0]]
				try:
					os.makedirs(outdir)
				except FileExistsError:
					cmd = ['WARNING: Output directory already exists;'] + cmd
				print(" ".join(cmd))
		elif i[1].endswith("gz"):
			print("unzip vcf file: ", i[1])		
	else:
		print(i[1], "does not exist") 