#!/usr/bin/env python
"""This script generates commands for mapping fastq files with 'speedseq' as well as subsequent commands for merging bams within a sample.  We map all fastqs corresponding to the same sample separately and then merge these after the fact in order to properly specify reads groups for duplicate removal.  These commands can be executed in parallel using GNU parallel or a task array.  Expects paired-end data for all samples and that such data is in separate files for forward and reverse reads.  For the mapping step, all bams within a sample will be written to a temporary folder named after the sample.  During the merging step, all bams will be merged within the temporary folder and written to the final output directory.  
 Takes the following arguments:
    -f :  REQUIRED: Full path to directory with input fastq files'
    -k :  REQUIRED: Space-delimited file to key that links fastq files to sample names. fastq file names can be partial, but must be complete enough to uniquely identify the fastq file. One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
    -r :  REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path'
    -o :  REQUIRED: Full path to output directory in which the MERGED bam files will be written
    -c :  Number_of_cores
    -m :  number of gigabytes of memory
    -s :  path to the speedseq directory
"""

# Import all necessary modules
import os
import argparse
import time

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for mapping fastq files with "speedseq" as well as subsequent commands for merging bams within a sample.  We map all fastqs corresponding to the same sample separately and then merge these after the fact in order to properly specify reads groups for duplicate removal.  These commands can be executed in parallel using GNU parallel or a task array.  Expects paired-end data for all samples and that such data is in separate files for forward and reverse reads.  For the mapping step, all bams within a sample will be written to a temporary folder named after the sample.  During the merging step, all bams will be merged within the temporary folder and written to the final output directory')
parser.add_argument('-f', type=str, metavar='fastq_directory', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-k', type=str, metavar='Sample_Fastq_Key', help='REQUIRED: Space-delimited file to key that links fastq files to sample names.  fastq file names can be partial, but must be complete enough to uniquely identify the fastq file.One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
parser.add_argument('-r', type=str, metavar='Reference_Path_Key', help='REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path')
parser.add_argument('-o', type=str, metavar='output_bam_directory', help='REQUIRED: Full path to output directory in which the MERGED bam files will be written')
parser.add_argument('-c', type=str, metavar='Number_of_mapping_cores', help="REQUIRED: Number_of_cores for mapping")
parser.add_argument('-m', type=str, metavar='Number_of_Gb', help="REQUIRED: number of gigabytes of memory for mapping")
parser.add_argument('-s', type=str, metavar='path_to_speedseq_directory', default="/home/hirschc1/pmonnaha/software/speedseq/")
args = parser.parse_args()


# Prepare necessary paths
if args.o.endswith("/") is False: args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if args.s.endswith("/") is False: args.s += "/"
    
# Read reference path information from Reference_Path_Key and store as dictionary where name of reference is the key and path to reference is the value
REFS = {}
with open(args.r, 'r') as ref_file:
    for line in ref_file:
        REFS[line.split()[0]] = line.split()[1]

# Read samples from the Sample_Fastq_Key and store in a dictionary where sample names are key and fastq names are values
samps = {}
with open(args.k, 'r') as file:
    for line in file:
        samps[line.split()[0]] = line.split()[1:]

# Main loop to write commands
for samp,fqs in samps.items():
    for ref, refpath in REFS.items():
        out_prefix = samp + "_" + ref
        BAM_string = ""
        for fq in fqs:
            lib = "b" # According to CNH, all fastqs within a sample that begin with numbered code are from same library.  Fastqs beginning with letter (corresponding to sample name) are sequences they got from elsewhere
            if fq[0].isdigit():
                lib = "a"

            FQs = [s for s in os.listdir(args.f) if fq in s and "sing" not in s] # find absolute path to R1 and R2 fastqs based on the fastq name while excluding the singles file that 'sickle' would have produces
            FQs.sort() # This should sort the fastqs so that read 1 file appears before read 2
            
            # Check to make sure that the search for fastq files corresponding the name specified by 'fq' uniquely identifies a single paire of forward and reverse read files
            assert len(FQs) >=2, "Did not find fastq files for %r" % fq  
            assert len(FQs) <=2, "Found multiple sets of fastq files for %r" % fq

            tmp_name = fq + "_" + ref
            tmpdir = args.o + tmp_name + "/"

            # Write speedseq command.  
            print('mkdir -p ' + tmpdir + ';' + args.s + 'bin/speedseq align -M ' + args.m + ' -t ' + args.c + r' -R "@RG\tID:' + fq + r'\tSM:' + samp + r'\tLB:' + samp + lib + r'" -T ' + tmpdir + 'tmp -o ' + tmpdir + tmp_name + " " + refpath + " " + args.f + FQs[0] + " " + args.f + FQs[1])





