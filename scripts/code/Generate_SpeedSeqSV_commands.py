#!/usr/bin/env python
"""This script generates commands for discovering structural variants with 'speedseq sv'.  Speedseq sv implements lumpy (for discovery) along with svtyper and CNVnator (for genotyping)
 Takes the following arguments:
    -b :  REQUIRED: Full path to directory with input bam files
    -r :  REQUIRED: Space-delimited file to key that links reference fasta path AND reference gff path to reference names. FORMAT: Reference_name Reference.fasta Reference_NonGenic.bed
    -o :  REQUIRED: Full path to output directory in which the per-individual vcf files will be written
    -c :  Number_of_cores
    -s :  path to the speedseq directory
"""

# Import all necessary modules
import os
import argparse

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for discovering structural variants with speedseq sv.  Speedseq sv implements lumpy (for discovery) along with svtyper and CNVnator (for genotyping)')
parser.add_argument('-b', type=str, metavar='bam_directory', default="/panfs/roc/scratch/pmonnaha/Maize/widiv_bams/merged/", help='Full path to directory with input merged bam files')
parser.add_argument('-r', type=str, metavar='Reference_Path_Key', default="/home/hirschc1/pmonnaha/JobScripts/accessory/reference_paths.txt", help='Space-delimited file to key that links reference fasta path AND reference non-genic bed file path to reference names. Format: Reference_name Reference.fasta Reference_NonGenic.bed')
parser.add_argument('-o', type=str, metavar='output_directory', default="/panfs/roc/scratch/pmonnaha/Maize/Lumpy/", help='Full path to output directory in which the per-individual vcf files will be written')
parser.add_argument('-c', type=str, metavar='Number_of_cores', default="8", help="Number_of_cores")
parser.add_argument('-s', type=str, metavar='path_to_speedseq_directory', default="/home/hirschc1/pmonnaha/software/speedseq/")
parser.add_argument('-w', type=str, metavar='window_size', default="300")
args = parser.parse_args()


# Prepare necessary paths
if args.o.endswith("/") is False: args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if args.s.endswith("/") is False: args.s += "/"
if args.b.endswith("/") is False: args.b += "/"
    
# Read reference path information from Reference_Path_Key and store as dictionary where name of reference is the key and path to reference is the value
REFS = {}
with open(args.r, 'r') as ref_file:
    for line in ref_file:
        REFS[line.split()[0]] = line.split()[1:]

# Main loop to write commands
for file in os.listdir(args.b):
    if file.endswith(".bam") and "disc" not in file and "splt" not in file:
        bam_name = file.strip(".bam")
        ref = bam_name.split("_")[1]

        # Print speedseq sv command.  
        cmd = args.s + "bin/speedseq sv -B " + args.b + file + " -S " + args.b + bam_name + ".splt.bam -D " + args.b + bam_name + ".disc.bam -R " + REFS[ref][0] + " -x " + REFS[ref][1] + " -o " + args.o + bam_name + " -K " + args.s + "bin/speedseq." + ref + ".config -t " + args.c + " -T " + args.o + bam_name + "_tmp -w " + args.w + " -v -d -P -g -k"

        if os.path.exists(args.o + bam_name + ".sv.vcf.gz"):
            cmd = "Error: Output file already exists; " + cmd

        print(cmd)


