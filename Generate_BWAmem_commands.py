import os, sys, argparse, subprocess, csv


parser = argparse.ArgumentParser(description='')

parser.add_argument('-fastqdir', type=str, metavar='fastq_ca', default='fastq_ca/', help='REQUIRED: Full path to directory with input fastq files')

parser.add_argument('-S', type=str, metavar='Sample_Fastq_Key', help='REQUIRED: File to key that links fastq files to sample names')
parser.add_argument('-o', type=str, metavar='aligned_dir', default='aligned/', help='Realtive path to output directory [aligned/]')

parser.add_argument('-oe', type=str, metavar='stdout_and_stderr_directory')
# parser.add_argument('-c', type=str, metavar='number_of_processors')

parser.add_argument('-s', type=str, metavar='samples', default='all')

parser.add_argument('-t', type=str, metavar='time', default='24:00:00')

parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')

args = parser.parse_args()

if args.o.endswith("/") is False:
	args.o += "/"

samps = {}
with open(args.S, 'r') as file:
	for line in file:
		if args.s != 'all':
			if line.split()[0] in args.s.split(","):
				samps[line.split()[0]] = line.split()[1:]
		else:
			samps[line.split()[0]] = line.split()[1:]

REFS = {'B73v4': '/home/hirschc1/pmonnaha/misc-files/B73_chr1-10.fasta', 'PH207': '/home/hirschc1/pmonnaha/misc-files/PH207_chr1-10.fasta', 'W22v12': '/home/hirschc1/pmonnaha/misc-files/W22_chr1-10.fasta','PHB47': '/home/hirschc1/pmonnaha/misc-files/PHB47_chr1-10.fasta'}
tmpdir = "/scratch.local/pmonnaha/"
# print(samps)
for samp,fqs in samps.items():
	for ref, refpath in REFS.items():
		out_prefix = samp + "_" + ref
		sh_file = open(args.o + out_prefix + '.sh', 'w')
		sh_file.write('#!/bin/bash -e\n'+
                  '#PBS -l mem=62gb,nodes=1:ppn=24,walltime=' + args.t + '\n'+
                  '#PBS -A hirschc1\n'+
                  '#PBS -q mesabi\n'+
                  '#PBS -o ' + args.oe + samp + "_" + ref + '.out\n'+ 
                  '#PBS -e ' + args.oe + samp + "_" + ref + '.err\n'+
                  '#PBS -M pmonnaha@umn.edu\n\n'+
                  'module load bwa\n'+
                  'module load samtools\n\n'+
                  'mkdir -p ' + tmpdir + '\n\n')
		FQ_string = ""
		for fq in fqs: # Need to see if the different fqs correspond to different libraries or if it was just the same library sequenced multiple times.
			lib = "b" # According to CNH, all fastqs within a sample that begin with numbered code are from same library.  Fastqs beginning with letter (corresponding to sample name) are sequences they got from elsewhere
			if fq[0].isdigit():
				lib = "a"
			FQs = [s for s in os.listdir(args.fastqdir) if fq in s and "singles" not in s] # find absolute path to R1 and R2 fastqs based on the fastq name
			# print(samp,fqs)
			# print(FQs)
			FQs.sort()
			if len(FQs) >= 2 and len(FQs) <= 3:
				fq_name = fq + "_" + ref
				sh_file.write(r'bwa mem -t 24 -R "@RG\tID:' + fq + r'\tSM:' + samp + r'\tLB:' + samp + lib + r'" ' + refpath + " " + args.fastqdir + FQs[0] + " " + args.fastqdir + FQs[1] + ' | samtools view -@ 24 -Sbh - | samtools sort -n -@ 24 -T ' + tmpdir + fq_name + '.tmp > ' + tmpdir + fq_name + '.bam\n\n')
				FQ_string += tmpdir + fq_name + '.bam '

			else:
				print("Did not find fastq files for", fq, "from",samp)

		sh_file.write('samtools merge -n ' + tmpdir + samp + "_" + ref + '.bam ' + FQ_string + '\n\n'
					  'samtools view -@ 24 -h ' + tmpdir + samp + "_" + ref + ".bam | /panfs/roc/groups/14/hirschc1/pmonnaha/software/speedseq//bin/samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 --splitterFile >(samtools view -Sb - > " + tmpdir + out_prefix + ".splitters.bam) --discordantFile >(samtools view -Sb - > " + tmpdir + out_prefix + ".disc.bam) -o >(samtools view -Sb - > " + tmpdir + out_prefix + ".bam)\n\n" +

					  "samtools view -@ 24 -h -b -u " + tmpdir + out_prefix + ".bam | samtools sort -@ 24 -o " + args.o + out_prefix + ".bam -T " + tmpdir + out_prefix + "\n\n" + 

					  "samtools view -h -b -u " + tmpdir + out_prefix + ".splitters.bam | samtools sort -@ 24 -o " + args.o + out_prefix + ".splitters.bam -T " + tmpdir + out_prefix+ "\n\n" + 

					  "samtools view -h -b -u " + tmpdir + out_prefix + ".disc.bam | samtools sort -@ 24 -o " + args.o + out_prefix + ".disc.bam -T " + tmpdir + out_prefix + '\n\n'+

					  "samtools index " + args.o + out_prefix + ".bam\n"+
					  "samtools index " + args.o + out_prefix + ".splitters.bam\n"+
					  "samtools index " + args.o + out_prefix + ".disc.bam\n")
		sh_file.close()
		if args.print == 'false':
			cmd = ('qsub ' + args.o + out_prefix + '.sh')
			p = subprocess.Popen(cmd, shell=True)
			sts = os.waitpid(p.pid, 0)[1]
		else:
			file = open(args.o + out_prefix + '.sh','r')
			data = file.read()
			print(data)
		os.remove(args.o + out_prefix + '.sh')
