#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script generates commands for merging a new set of (unmerged) bams with a set of previously merged bams, using sambamba.
 Takes the following arguments:
    -u :  REQUIRED. unmerged_bam_directory. Full path to directory with input unmerged bam files
    -m :  REQUIRED. merged_bam_directory.  Full path to directory with input merged bam files
    -k :  REQUIRED. merge_Key.  Space-delimited file to key that contains the sample name (as named in the merged_bam_directory) and the bam prefixes of unmerged bams in subsequent columns.
    -o :  REQUIRED: Full path to output directory in which the MERGED bam files will be written
    -r :  REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path'
    -c :  Number_of_cores
    -s :  path to the speedseq directory
"""

# Import all necessary modules
import os
import argparse
import time
import pdb
import sys
from pathlib import Path
import shutil

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for merging a new set of (unmerged) bams with a set of previously merged bams, using sambamba.')
parser.add_argument('-u', type=str, metavar='unmerged_bam_directory', required=True, help='Full path to directory with input unmerged bam files')
parser.add_argument('-m', type=str, metavar='merged_bam_directory', required=True, help='Full path to directory with input merged bam files')
parser.add_argument('-k', type=str, metavar='merge_Key', required=True, help='Space-delimited file to key that contains the sample name (as named in the merged_bam_directory) and the bam prefixes of unmerged bams in subsequent columns.')
parser.add_argument('-o', type=str, metavar='output_bam_directory', required=True, help='Full path to output directory in which the new MERGED bam files will be written')
parser.add_argument('-r', type=str, metavar='Reference_Path_Key', required=True, help='Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path')
parser.add_argument('-c', type=str, metavar='Number_of_cores', default="1")
parser.add_argument('-s', type=str, metavar='path_to_speedseq_directory', default="~/speedseq/")
parser.add_argument('-d', type=str, metavar='delete_input', default='false', help='Do you want to delete the input files upon successful completion of the merging?')
args = parser.parse_args()

#check python version
assert sys.version_info[0] != 2, print("This script requires python3.  Do 'module load python/3.6.3'")

# Prepare necessary paths
if args.o.endswith("/") is False: args.o += "/"
if os.path.exists(args.o) is False: os.mkdir(args.o)
if args.s.endswith("/") is False: args.s += "/"

#check that sambamba executable is exists
home = str(Path.home())
sambamba_path = args.s.replace("~", home) + "/bin/sambamba"

assert shutil.which(sambamba_path), print("Did not find the sambamba executable at", sambamba_path, ".  Be sure to provide exact path to the speedseq repo that you cloned (see README).")

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
        missing = []
        found = 0
        cmd = ""
        if len(bams) > 0 and any(out_prefix in k for k in os.listdir(args.m)): # At least one bam is expected for this sample and the pre-existing merged bam (which should have the out_prefix as prefix) is in the merged bam directory
            for bam in bams: # loop over all component bams for same sample
                tmp_name = bam + "_" + ref
                # pdb.set_trace()
                if os.path.exists(args.u + tmp_name + "/" + tmp_name + ".bam"):
                    BAM_string += args.u + tmp_name + "/" + tmp_name + ".bam "
                    found += 1
                else:
                    missing.append(bam)

            mBAM_string = f"{args.m}/{out_prefix}.bam"

            fin_name = args.o + samp + "_" + ref + '.bam' # Name of final merged bam file

            # Write commands for merging bams within a sample.  Three commands per fastq are needed as each fastq produces a full, discordant-read, and split-read bam.
            cmd = f"{args.s}/bin/sambamba merge -t {args.c} {fin_name} {BAM_string} {mBAM_string} && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t {args.c} {fin_name} && {args.s}/bin/sambamba merge -t {args.c} {fin_name.replace('.bam','.disc.bam ')} {BAM_string.replace('.bam ', '.discordants.bam ')} {mBAM_string.replace('.bam', '.disc.bam')} && {args.s}/bin/sambamba index -t {args.c} {fin_name.replace('.bam', '.disc.bam')} && {args.s}/bin/sambamba merge -t {args.c} {fin_name.replace('.bam','.splt.bam ')} {BAM_string.replace('.bam ', '.splitters.bam ')} {mBAM_string.replace('.bam', '.splt.bam')} && {args.s}/bin/sambamba index -t {args.c} {fin_name.replace('.bam', '.splt.bam')}"
            if args.d == "true": # Delete input files?
                cmd += " && rm " + BAM_string + " && rm " + BAM_string.replace(".bam", ".discordants.bam") + " && rm " + BAM_string.replace(".bam", ".splitters.bam") + " && rm " + mBAM_string.replace(".bam", ".disc.bam") + " && rm " + mBAM_string.replace(".bam", ".splt.bam")

            if missing != []:
                cmd = "Error: Did not find bams " + str(missing) + " for " + samp + ". " + cmd
            if os.path.exists(fin_name):
                cmd = "Error: " + fin_name + " already exists; " + cmd
            print(cmd)
        else:
            print(out_prefix, "not found")

