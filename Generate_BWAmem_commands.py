import os, sys, argparse, subprocess, csv, time


parser = argparse.ArgumentParser(description='')

parser.add_argument('-fastqdir', type=str, metavar='fastq_ca', default='fastq_ca/', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-S', type=str, metavar='Sample_Fastq_Key', help='REQUIRED: File to key that links fastq files to sample names')
parser.add_argument('-o', type=str, metavar='aligned_dir', default='aligned/', help='Realtive path to output directory [aligned/]')
parser.add_argument('-c', type=str, metavar='command_directory', help="directory in which to put the three command files that are produced")
parser.add_argument('-C', type=str, metavar='Number_of_cores', help="Number_of_cores for mapping")

parser.add_argument('-s', type=str, metavar='samples', default='all')
parser.add_argument('-n', type=str, metavar='run_name', default='-99')

args = parser.parse_args()


if args.o.endswith("/") is False:
	args.o += "/"
if os.path.exists(args.o + "Unmerged") is False:
	os.mkdir(args.o + "Unmerged/")
tmpdir = args.o + "Unmerged/"

if args.c.endswith("/") is False:
	args.c += "/"

ID = str(time.time())
if args.n != "-99":
	ID = args.n

bwa_out = open(args.c + "bwa_mem_commands_" + ID + '.txt', 'w')
ms_out = open(args.c + "merge_and_split_commands_" + ID + '.txt', 'w')
s_out = open(args.c + "sort_commands_" + ID + '.txt', 'w')

REFS = {'B73v4': '/home/hirschc1/pmonnaha/misc-files/B73_chr1-10.fasta', 'PH207': '/home/hirschc1/pmonnaha/misc-files/PH207_chr1-10.fasta', 'W22v12': '/home/hirschc1/pmonnaha/misc-files/W22_chr1-10.fasta','PHB47': '/home/hirschc1/pmonnaha/misc-files/PHB47_chr1-10.fasta'}


samps = {}
with open(args.S, 'r') as file:
	for line in file:
		if args.s != 'all':
			if line.split()[0] in args.s.split(","):
				samps[line.split()[0]] = line.split()[1:]
		else:
			samps[line.split()[0]] = line.split()[1:]


for samp,fqs in samps.items():
	for ref, refpath in REFS.items():
		out_prefix = samp + "_" + ref
		FQ_string = ""
		for fq in fqs: # Need to see if the different fqs correspond to different libraries or if it was just the same library sequenced multiple times.
			lib = "b" # According to CNH, all fastqs within a sample that begin with numbered code are from same library.  Fastqs beginning with letter (corresponding to sample name) are sequences they got from elsewhere
			if fq[0].isdigit():
				lib = "a"
			FQs = [s for s in os.listdir(args.fastqdir) if fq in s and "sing" not in s] # find absolute path to R1 and R2 fastqs based on the fastq name
			sings = [s for s in os.listdir(args.fastqdir) if fq in s and "sing" in s]

			FQs.sort()
			if len(FQs) >= 2 and len(FQs) <= 3:
				fq_name = fq + "_" + ref
				bwa_out.write(r'bwa mem -M -t ' + str(args.C) + r' -R "@RG\tID:' + fq + r'\tSM:' + samp + r'\tLB:' + samp + lib + r'" ' + refpath + " " + args.fastqdir + FQs[0] + " " + args.fastqdir + FQs[1] + ' -U ' + args.fastqdir + sings[0] + ' | samtools view -@ ' + str(args.C) + ' -Sbh - | samtools sort -n -@ ' + str(args.C) + ' -T ' + tmpdir + fq_name + '.tmp > ' + tmpdir + fq_name + '.bam\n')

				FQ_string += tmpdir + fq_name + '.bam '

			else:
				print("Did not find fastq files for", fq, "from",samp)

		ms_out.write('samtools merge -n ' + args.o + samp + "_" + ref + '.bam ' + FQ_string + ' && rm ' + FQ_string + 
					  '; samtools view -h ' + args.o + samp + "_" + ref + ".bam | /panfs/roc/groups/14/hirschc1/pmonnaha/software/speedseq//bin/samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 --splitterFile >(samtools view -Sb - | samtools sort -o " + args.o + out_prefix + ".splitters.bam -T " + args.o + out_prefix + ".split) --discordantFile >(samtools view -Sb - | samtools sort -o " + args.o + out_prefix + ".disc) -o >(samtools view -Sb - > " + args.o + out_prefix + ".unsort.bam)\n")


		s_out.write("samtools view -@ 24 -h -b -u " + args.o + out_prefix + ".unsort.bam | samtools sort -@ 24 -o " + args.o + out_prefix + ".bam -T " + args.o + out_prefix + " && rm " + args.o + out_prefix + ".unsort.bam; samtools index " + args.o + out_prefix + ".bam\n")

