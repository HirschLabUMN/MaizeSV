#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script generates commands for merging files with sambamba.  These commands can be executed in parallel using GNU parallel or a task array.  Be sure to use the same Sample_Fastq_Key (-s) and output_directory (-o) as was used with Generate_SpeedSeq_commands.py
 Takes the following arguments:
    -s :  REQUIRED: Space-delimited file to key that links fastq files to sample names. fastq file names can be partial, but must be complete enough to uniquely identify the fastq file. One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
    -o :  REQUIRED: Full path to output directory in which the MERGED bam files will be written
    -v :  REQUIRED: merged vcf file resulting from svtools lmerge.  The portion of the prefix denoting the reference should match the specification in the bam file names.
    -w :  window size in which CNVnator will calculate depth
"""

# Import all necessary modules
import os
import argparse
import time

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for merging files with sambamba.  These commands can be executed in parallel using GNU parallel or a task array.  Be sure to use the same Sample_Fastq_Key (-s) and output_directory (-o) as was used with Generate_SpeedSeq_commands.py')
parser.add_argument('-b', type=str, metavar='bam_directory', required=True, help='Full path to directory with input unmerged bam files')
parser.add_argument('-o', type=str, metavar='output_directory', help='REQUIRED: Full path to output directory in which the gt and cn subdirectories are located')
parser.add_argument('-v', type=str, metavar='merged_vcf', required=True, help='merged vcf file resulting from svtools lmerge')
parser.add_argument('-c', type=str, metavar='coord_file', required=True, help='path to coordinates file created via "create_coordinates -i merged.vcf -o coord_file"')
parser.add_argument('-w', type=str, metavar='window_size', default="100", help='window size in which CNVnator will calculate depth')
args = parser.parse_args()


# Prepare necessary paths
if args.o.endswith("/") is False: args.o += "/"
if os.path.exists(args.o) is False: print("Output directory does not exist.  Use same output directory that was used for Speedseq")
if args.s.endswith("/") is False: args.s += "/"
if args.b.endswith("/") is False: args.b += "/"

# Main loop to write commands
for bam in os.listdir(args.b):
    if bam.endswith(".bam") and "splt" not in bam and "disc" not in bam:
        ref = bam.split("_")[1]
        sample = bam.split("_")[0]
        if ref in args.v:
            print("zcat " + args.v + " | vawk --header '{  $6=\".\"; print }' | svtools genotype -B " + bam + " -l " + bam + r".json | sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - > " + args.o + "gt/" + bam.strip(".bam") + ".vcf && svtools copynumber --cnvnator cnvnator-multi -s " + sample + " -w " + args.w + " -r " + args.b + bam.strip(".bam") + " -c " + args.c + " -i " + args.o + "gt/" + bam.strip(".bam") + ".vcf" > args.o + "cn/" + bam.strip(".bam") + ".vcf")
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

mkdir -p gt

zcat merged.vcf.gz \
| vawk --header '{  $6="."; print }' \
| svtools genotype \
  -B NA12877.bam \
  -l NA12877.bam.json \
| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
> gt/NA12877.vcf

mkdir -p cn

svtools copynumber \
--cnvnator cnvnator-multi \
-s NA12877 \
-w 100 \
-r /temp/cnvnator-temp/NA12877.bam.hist.root \
 -c coordinates \
 -i gt/NA12877.vcf \
> cn/NA12877.vcf


        # Write commands for merging bams within a sample.  Three commands per fastq are needed as each fastq produces a full, discordant-read, and split-read bam.

