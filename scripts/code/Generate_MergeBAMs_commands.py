#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script generates commands for merging files with sambamba.  These commands can be executed in parallel using GNU parallel or a task array.  Be sure to use the same Sample_Fastq_Key (-s) and output_directory (-o) as was used with Generate_SpeedSeq_commands.py
 Takes the following arguments:
    -s :  REQUIRED: Space-delimited file to key that links fastq files to sample names. fastq file names can be partial, but must be complete enough to uniquely identify the fastq file. One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
    -o :  REQUIRED: Full path to output directory in which the MERGED bam files will be written
    -r :  REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path'
    -c :  Number_of_cores
    -s :  path to the speedseq directory
"""

# Import all necessary modules
import os
import argparse
import time

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for merging files with sambamba.  These commands can be executed in parallel using GNU parallel or a task array.  Be sure to use the same Sample_Fastq_Key (-s) and output_directory (-o) as was used with Generate_SpeedSeq_commands.py')
parser.add_argument('-b', type=str, metavar='bam_directory', required=True, help='Full path to directory with input unmerged bam files')
parser.add_argument('-k', type=str, metavar='Sample_Fastq_Key', help='REQUIRED: Space-delimited file to key that links fastq files to sample names.  fastq file names can be partial, but must be complete enough to uniquely identify the fastq file.One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
parser.add_argument('-o', type=str, metavar='output_bam_directory', help='REQUIRED: Full path to output directory in which the MERGED bam files will be written')
parser.add_argument('-c', type=str, metavar='Number_of_cores', default="1")
parser.add_argument('-s', type=str, metavar='path_to_speedseq_directory', default="/home/hirschc1/pmonnaha/software/speedseq/")
parser.add_argument('-r', type=str, metavar='Reference_Path_Key', help='REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path')
args = parser.parse_args()


# Prepare necessary paths
if args.o.endswith("/") is False: args.o += "/"
if os.path.exists(args.o) is False: print("Output directory does not exist.  Use same output directory that was used for Speedseq")
if args.s.endswith("/") is False: args.s += "/"

# Read samples from the Sample_Fastq_Key and store in a dictionary where sample names are key and fastq names are values
samps = {}
with open(args.k, 'r') as file:
    for line in file:
        samps[line.split()[0]] = line.split()[1:]

# Read reference path information from Reference_Path_Key and store as dictionary where name of reference is the key and path to reference is the value
REFS = {}
with open(args.r, 'r') as ref_file:
    for line in ref_file:
        REFS[line.split()[0]] = line.split()[1]

# Main loop to write commands
for samp,bams in samps.items():
    for ref, refpath in REFS.items():
        out_prefix = samp + "_" + ref
        BAM_string = ""
        for bam in bams:
            tmp_name = bam + "_" + ref
            tmpdir = args.o + tmp_name + "/"
            bam_name = [s for s in os.listdir(args.b) if bam in s and "splt" not in s and "disc" not in s]
            
            assert len(bam_name) == 1, "Found multiple files corresponding to %r" % bam

            BAM_string += tmpdir + bam_name

            # Make sure that all of the expected (full, discordant, and split reads) bams exist
            error = True
            if os.path.exists(tmpdir + tmp_name + '.bam') and os.path.exists(tmpdir + tmp_name + '.discordants.bam') and os.path.exists(tmpdir + tmp_name + '.splitters.bam'):
                error = False
            assert error is False, "Did not find all bams (full, discordant, and split reads) for fastq file: %r" % bam

        fin_name = args.o + samp + "_" + ref + '.bam' # Name of final merged bam file

        # Write commands for merging bams within a sample.  Three commands per fastq are needed as each fastq produces a full, discordant-read, and split-read bam.
        print(args.s + 'bin/sambamba merge -t ' + args.c + " " + fin_name + ' ' + BAM_string + ' && rm ' + BAM_string + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + args.c + " " + fin_name)
        print(args.s + 'bin/sambamba merge -t ' + args.c + " " + fin_name.replace('.bam','.disc.bam ') + BAM_string.replace(".bam", ".discordants.bam") + ' && rm ' + BAM_string.replace(".bam", ".discordants.bam") + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + args.c + " " + fin_name.replace('.bam', '.disc.bam'))
        print(args.s + 'bin/sambamba merge -t ' + args.c + " " + fin_name.replace('.bam','.splt.bam ') + BAM_string.replace(".bam", ".splitters.bam") + ' && rm ' + BAM_string.replace(".bam", ".splitters.bam") + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + args.c + " " + fin_name.replace('.bam', '.splt.bam'))

