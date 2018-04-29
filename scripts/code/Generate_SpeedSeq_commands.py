#!/usr/bin/env python
"""This script generates commands for mapping fastq files with 'speedseq' as well as subsequent commands for merging bams within a sample.  We map all fastqs corresponding to the same sample separately and then merge these after the fact in order to properly specify reads groups for duplicate removal.  These commands can be executed in parallel using GNU parallel or a task array.  Expects paired-end data for all samples and that such data is in separate files for forward and reverse reads.  For the mapping step, all bams within a sample will be written to a temporary folder named after the sample.  During the merging step, all bams will be merged within the temporary folder and written to the final output directory.  
 Takes four arguments:
    -f :  REQUIRED: Full path to directory with input fastq files'
    -s :  REQUIRED: Space-delimited file to key that links fastq files to sample names. fastq file names can be partial, but must be complete enough to uniquely identify the fastq file. One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
    -r :  REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path'
    -o :  REQUIRED: Full path to output directory in which the MERGED bam files will be written
    -c :  REQUIRED: directory in which to put the two command files that are produced
    -C :  REQUIRED: Number_of_cores for mapping
    -M :  REQUIRED: number of gigabytes of memory for mapping
    -C2:  Number_of_cores for merging.  If not specified, will use value for -C
    -n :  Name of the mapping/merging run.  Added to the command files for personal reference
"""

# Import all necessary modules
import os
import argparse
import time

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser(description='This script generates commands for mapping fastq files with "speedseq" as well as subsequent commands for merging bams within a sample.  We map all fastqs corresponding to the same sample separately and then merge these after the fact in order to properly specify reads groups for duplicate removal.  These commands can be executed in parallel using GNU parallel or a task array.  Expects paired-end data for all samples and that such data is in separate files for forward and reverse reads.  For the mapping step, all bams within a sample will be written to a temporary folder named after the sample.  During the merging step, all bams will be merged within the temporary folder and written to the final output directory')
parser.add_argument('-f', type=str, metavar='fastq_directory', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-s', type=str, metavar='Sample_Fastq_Key', help='REQUIRED: Space-delimited file to key that links fastq files to sample names.  fastq file names can be partial, but must be complete enough to uniquely identify the fastq file.One sample per line. Format: Sample_name fastq1_prefix fastq2_prefix fastq3_prefix')
parser.add_argument('-r', type=str, metavar='Reference_Path_Key', help='REQUIRED: Space-delimited file to key that links reference fasta path to reference names to be used in BAM naming. Format: Reference_name Reference_path')
parser.add_argument('-o', type=str, metavar='output_bam_directory', help='REQUIRED: Full path to output directory in which the MERGED bam files will be written')
parser.add_argument('-c', type=str, metavar='command_directory', default='/home/hirschc1/pmonnaha/JobScripts/accessory/', help="REQUIRED: directory in which to put the two command files that are produced")
parser.add_argument('-C', type=str, metavar='Number_of_mapping_cores', help="REQUIRED: Number_of_cores for mapping")
parser.add_argument('-M', type=str, metavar='Number_of_Gb', help="REQUIRED: number of gigabytes of memory for mapping")
parser.add_argument('-C2', type=str, metavar='Number_of_merging_cores', default="-99", help="Number_of_cores for merging; If not specified, will use value for -C")
parser.add_argument('-n', type=str, metavar='run_name', default='-99')
args = parser.parse_args()


# Prepare necessary directories for output
if args.o.endswith("/") is False: args.o += "/"
if args.c.endswith("/") is False: args.c += "/" 
if os.path.exists(args.o) is False: os.mkdir(args.o)

# If a name for the command files is not provided, use the time as a unique identifier
ID = str(time.time())
if args.n != "-99": ID = args.n

# Set number of cores to be used for merging:
if args.C2 == "-99":
    C2 = args.C
else:
    C2 = args.C2
    
# Open files that will store the commands
ss_out = open(args.c + "speedseq_commands_" + ID + '.txt', 'w')
ms_out = open(args.c + "merge_and_split_commands_" + ID + '.txt', 'w')
s_out = open(args.c + "sort_commands_" + ID + '.txt', 'w')

# Read reference path information from Reference_Path_Key and store as dictionary where name of reference is the key and path to reference is the value
REFS = {}
with open(args.r, 'r') as ref_file:
    for line in ref_file:
        REFS[line.split()[0]] = line.split()[1]

# Read samples from the Sample_Fastq_Key and store in a dictionary where sample names are key and fastq names are values
samps = {}
with open(args.s, 'r') as file:
    for line in file:
        samps[line.split()[0]] = line.split()[1:]

# Main loop to write commands
unfound = [] # Store fastq file names specified in Sample_fastq_key, but not found in fastq directory
bad_names = [] # Store fastq file names specified in Sample_fastq_key that do not uniquely identify a single set of forward and reverse reads
for samp,fqs in samps.items():
    for ref, refpath in REFS.items():
        out_prefix = samp + "_" + ref
        BAM_string = ""
        for fq in fqs: # Need to see if the different fqs correspond to different libraries or if it was just the same library sequenced multiple times.
            lib = "b" # According to CNH, all fastqs within a sample that begin with numbered code are from same library.  Fastqs beginning with letter (corresponding to sample name) are sequences they got from elsewhere
            if fq[0].isdigit():
                lib = "a"

            FQs = [s for s in os.listdir(args.f) if fq in s and "sing" not in s] # find absolute path to R1 and R2 fastqs based on the fastq name while excluding the singles file that 'sickle' would have produces
            FQs.sort() # This should sort the fastqs so that read 1 file appears before read 2
            if len(FQs) == 2: # Only write command if 
                tmp_name = fq + "_" + ref
                tmpdir = args.o + tmp_name + "/"

                # Write speedseq command.  
                ss_out.write('mkdir -p ' + tmpdir + r';/home/hirschc1/pmonnaha/software/speedseq/bin/speedseq align -M ' + args.M + ' -t ' + args.C + r' -R "@RG\tID:' + fq + r'\tSM:' + samp + r'\tLB:' + samp + lib + r'" -T ' + tmpdir + 'tmp -o ' + tmpdir + tmp_name + " " + refpath + " " + args.f + FQs[0] + " " + args.f + FQs[1] + "\n")

                BAM_string += tmpdir + tmp_name + '.bam '

            elif len(FQs) >2:
                bad_names.append("Found multiple sets of fastq files for " + samp + " corresponding to fastq file name: " + fq)
            else:
                unfound.append("Did not find fastq files for " + fq + " from " + samp)

        fin_name = args.o + samp + "_" + ref + '.bam' # Name of final merged bam file

        # Write commands for merging bams within a sample.  Three commands per fastq are needed as each fastq produces a full, discordant-read, and split-read bam.
        ms_out.write('/home/hirschc1/pmonnaha/software/speedseq/bin/sambamba merge -t ' + C2 + " " + fin_name + ' ' + BAM_string + ' && rm ' + BAM_string + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + C2 + " " + fin_name)
        ms_out.write('/home/hirschc1/pmonnaha/software/speedseq/bin/sambamba merge -t ' + C2 + " " + fin_name.replace('.bam','.disc.bam ') + BAM_string.replace(".bam", ".discordants.bam") + ' && rm ' + BAM_string.replace(".bam", ".discordants.bam") + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + args.C2 + " " + fin_name.replace('.bam', '.disc.bam'))
        ms_out.write('/home/hirschc1/pmonnaha/software/speedseq/bin/sambamba merge -t ' + C2 + " " + fin_name.replace('.bam','.splt.bam ') + BAM_string.replace(".bam", ".splitters.bam") + ' && rm ' + BAM_string.replace(".bam", ".splitters.bam") + " && /home/hirschc1/pmonnaha/software/speedseq/bin/sambamba index -t " + C2 + " " + fin_name.replace('.bam', '.splt.bam'))

# Print out the fastq file names that we were unable to find (missing_fastq) or fastq names that did not uniquely identify a set of forward and reverse reads (bad_names)
for missing_fastq in unfound: print(missing_fastq)
for bad_name in bad_names: print(bad_name)


