#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script generates commands for svtools genotype and svtools copynumber.  These commands are meant to be run subsequent to generating a merged vcf via svtools lsort followed by lmerge.  copynumber must be run subsequently to genotype and will run only if the latter executes successfully.  The corresponding shell script to run these commands is SVtools_Genotype.sh.
 Takes the following arguments:
    -b :  REQUIRED: Full path to directory with input bam files used for speedseq sv
    -o :  REQUIRED: Full path to output directory.  Results for svtools genotype and svtools copynumber will be stored within subdirectories gt and cn, respectively.  If these already exist and contain output from a previous run, these files will be overwritten.
    -v :  REQUIRED: merged vcf file resulting from svtools lmerge.  The portion of the prefix denoting the reference should match the specification in the bam file names.
    -c :  REQUIRED: path to coordinates file created via "create_coordinates -i merged.vcf -o coord_file"
    -w :  window size in which CNVnator will calculate depth
    -s :  path to speedseq directory
"""

# Import all necessary modules
import os
import argparse
import time

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for svtools genotype and svtools copynumber.  These commands are meant to be run subsequent to generating a merged vcf via svtools lsort followed by lmerge.  copynumber must be run subsequently to genotype and will run only if the latter executes successfully.  The corresponding shell script to run these commands is SVtools_Genotype.sh.')
parser.add_argument('-b', type=str, metavar='bam_directory', required=True, help='Full path to directory with input bam files used for speedseq sv')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Full path to output directory.  Results for svtools genotype and svtools copynumber will be stored within subdirectories gt and cn, respectively.  If these already exist and contain output from a previous run, these files will be overwritten')
# parser.add_argument('-l', type=str, metavar='lumpy_output_directory', required=True, help='Directory that contians the output from lumpy/speedseq sv')
parser.add_argument('-v', type=str, metavar='merged_vcf', required=True, help='merged vcf file resulting from svtools lmerge.  The portion of the prefix denoting the reference should match the specification in the bam file names.')
parser.add_argument('-c', type=str, metavar='coord_file', required=True, help='path to coordinates file created via "create_coordinates -i merged.vcf -o coord_file"')
parser.add_argument('-w', type=str, metavar='window_size', default="100", help='window size in which CNVnator will calculate depth')
parser.add_argument('-s', type=str, metavar='path_to_speedseq_directory', default="/home/hirschc1/pmonnaha/software/speedseq/")
parser.add_argument('-t', type=str, metavar='number_of_threads', default="1")
args = parser.parse_args()

vcf = args.v

# Prepare necessary paths
if args.o.endswith("/") is False: args.o += "/"
if os.path.exists(args.o) is False: print("Output directory does not exist")
if os.path.exists(args.o + "gt") or os.path.exists(args.o + "cn"): # Check if necessary subdirectory structure already exists
    print("Output subdirectories already exist.  !!Running the resulting commands will overwrite files present in these subdirectories!!")
else:
    os.mkdir(args.o + "cn")
if args.s.endswith("/") is False: args.s += "/"
if args.b.endswith("/") is False: args.b += "/"

# Main loop to write commands
for bam in os.listdir(args.b):
    if bam.endswith(".bam") and "splt" not in bam and "disc" not in bam: # We only want the full bams, not the discordant or split read bams
        ref = bam.split("_")[1].split(".")[0] 
        sample = bam.split("_")[0]
        if ref in vcf: # Only use the bams that pertain to the same reference genome in the merged vcf file
            print("/home/hirschc1/pmonnaha/software/speedseq/bin/cnvnator_wrapper.py -t " + str(args.t) + " -w " + args.w + " -b " + args.b + bam + " -c " + "/home/hirschc1/pmonnaha/references/" + ref + "_chroms/ -T " + args.o + bam.strip(".bam") + "_tmp/cnvnator-temp -o " + args.o + bam.strip(".bam") + "_tmp/cnvnator-temp/" + bam.strip(".bam") + "_" + args.w + ".readdepth --cnvnator " + args.s + "bin/cnvnator --samtools " + args.s + "src/samtools-1.3.1/samtools && svtools copynumber --cnvnator " + args.s + "bin/cnvnator -s " + sample + " -w " + args.w + " -r " + args.o + bam.strip(".bam") + "_tmp/cnvnator-temp/" + bam + ".hist.root -c " + args.c + " -i " + args.o + "gt/" + bam.strip(".bam") + ".vcf > " + args.o + "cn/" + bam.strip(".bam") + ".vcf")

