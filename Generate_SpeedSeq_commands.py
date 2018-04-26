#!/usr/bin/env python
"""This script generates commands to be used with for trimming adapters using cutadapt.  The resulting commands can be executed in parallel using GNU parallel or a task array.  Should only be used for paired-end data.  Output is meant to be .
 Takes four arguments:
    -r  :  REQUIRED. Full path to file with a list of read paths, one per line.  IMPORTANT: if paired reads are in separate files (i.e. not interleaved), create two files: one for forward reads and one for reverse reads. Pass the file containing paths to reverse read files to -r2
    -r2 : file containing paths to reverse read files
    -o  : Directory where the output of cutadapt will be written
    -c  : path to the cutadapt executable.  The most recent version (1.16) is required, which is not available on the cluster.  Download and install 1.16 and provide the path to cutadapt using this flag.
    -s  : path to sickle executable
"""

import os
import argparse
import time

parser = argparse.ArgumentParser(description='')

parser.add_argument('-fastqdir', type=str, metavar='fastq_ca', default='fastq_ca/', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-S', type=str, metavar='Sample_Fastq_Key', help='REQUIRED: Space-delimited file to key that links fastq files to sample names. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
parser.add_argument('-o', type=str, metavar='aligned_dir', default='aligned/', help='Realtive path to output directory [aligned/]')
parser.add_argument('-c', type=str, metavar='command_directory', help="directory in which to put the three command files that are produced")
parser.add_argument('-C', type=str, metavar='Number_of_mapping_cores', help="Number_of_cores for mapping")
parser.add_argument('-C2', type=str, metavar='Number_of_merging_cores', help="Number_of_cores for merging")
parser.add_argument('-M', type=str, metavar='Number_of_Gb', help="number of gigabytes of memory")

parser.add_argument('-s', type=str, metavar='samples', default='all')
parser.add_argument('-n', type=str, metavar='run_name', default='-99')

args = parser.parse_args()


if args.o.endswith("/") is False:
    args.o += "/"

if os.path.exists(args.o + "Unmerged") is False:
    os.mkdir(args.o + "Unmerged/")

if args.c.endswith("/") is False:
    args.c += "/"

ID = str(time.time())
if args.n != "-99":
    ID = args.n

ss_out = open(args.c + "speedseq_commands_" + ID + '.txt', 'w')
ms_out = open(args.c + "merge_and_split_commands_" + ID + '.txt', 'w')
s_out = open(args.c + "sort_commands_" + ID + '.txt', 'w')

REFS = {'B73v4':'/home/hirschc1/pmonnaha/misc-files/B73_chr1-10.fasta', 'PH207':'/home/hirschc1/pmonnaha/misc-files/PH207_chr1-10.fasta', 'W22v12':'/home/hirschc1/pmonnaha/misc-files/W22_chr1-10.fasta','PHB47':'/home/hirschc1/pmonnaha/misc-files/PHB47_chr1-10.fasta'}

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
        BAM_string = ""
        for fq in fqs: # Need to see if the different fqs correspond to different libraries or if it was just the same library sequenced multiple times.
            lib = "b" # According to CNH, all fastqs within a sample that begin with numbered code are from same library.  Fastqs beginning with letter (corresponding to sample name) are sequences they got from elsewhere
            if fq[0].isdigit():
                lib = "a"
            FQs = [s for s in os.listdir(args.fastqdir) if fq in s and "sing" not in s] # find absolute path to R1 and R2 fastqs based on the fastq name

            FQs.sort()
            if len(FQs) >= 2 and len(FQs) <= 3:
                fq_name = fq + "_" + ref
                outdir = args.o + "Unmerged/"
                tmpdir = args.o + "Unmerged/" + fq_name + "/"
                ss_out.write('mkdir -p ' + tmpdir + r';/home/hirschc1/pmonnaha/software/speedseq/bin/speedseq align -M ' + str(args.M) + ' -t ' + str(args.C) + r' -R "@RG\tID:' + fq + r'\tSM:' + samp + r'\tLB:' + samp + lib + r'" -T ' + tmpdir + ' -o ' + tmpdir + fq_name + " " + refpath + " " + args.fastqdir + FQs[0] + " " + args.fastqdir + FQs[1] + "\n")

                BAM_string += tmpdir + fq_name + '.bam '

            else:
                print("Did not find fastq files for", fq, "from",samp)

        fin_name = args.o + samp + "_" + ref + '.bam'
        ms_out.write('/home/hirschc1/pmonnaha/software/speedseq/bin/sambamba merge -t ' + str(args.C2) + " " + fin_name + ' ' + BAM_string + ' && rm ' + BAM_string + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + str(args.C2) + " " + fin_name)
        ms_out.write('/home/hirschc1/pmonnaha/software/speedseq/bin/sambamba merge -t ' + str(args.C2) + " " + fin_name.replace('.bam','.disc.bam ') + BAM_string.replace(".bam", ".discordants.bam") + ' && rm ' + BAM_string.replace(".bam", ".discordants.bam") + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + str(args.C2) + " " + fin_name.replace('.bam', '.disc.bam'))
        ms_out.write('/home/hirschc1/pmonnaha/software/speedseq/bin/sambamba merge -t ' + str(args.C2) + " " + fin_name.replace('.bam','.splt.bam ') + BAM_string.replace(".bam", ".splitters.bam") + ' && rm ' + BAM_string.replace(".bam", ".splitters.bam") + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + str(args.C2) + " " + fin_name.replace('.bam', '.splt.bam'))


        s_out.write("samtools view -@ 24 -h -b -u " + args.o + out_prefix + ".unsort.bam | samtools sort -@ 24 -o " + args.o + out_prefix + ".bam -T " + args.o + out_prefix + " && rm " + args.o + out_prefix + ".unsort.bam; samtools index " + args.o + out_prefix + ".bam\n")

