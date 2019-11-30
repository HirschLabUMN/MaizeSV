#!/usr/bin/env python
"""This script generates commands to be used with for trimming adapters using cutadapt.  The resulting commands can be executed in parallel using GNU parallel or a task array.  Should only be used for paired-end data.  Output is meant to be .
 Takes four arguments:
    -r  :  REQUIRED. Full path to file with a list of read paths, one per line.  IMPORTANT: if paired reads are in separate files (i.e. not interleaved), create two files: one for forward reads and one for reverse reads. Pass the file containing paths to reverse read files to -r2
    -r2 : file containing paths to reverse read files
    -o  : Directory where the output of cutadapt will be written
    -c  : path to the cutadapt executable.  The most recent version (1.16) is required, which is not available on the cluster.  Download and install 1.16 and provide the path to cutadapt using this flag.
    -s  : path to sickle executable
"""

# Import neeeded modules
import sys
import os
import pipes
import re
import subprocess
import argparse
from pathlib import Path
import shutil
import sys


def get_read_paths(read_path_file, read_path_file2):
    """Returns a list(s) which stores the path to each fastq file or pair of files.  These lists will be looped over and passed to build_command"""
    
    read_paths = []
    read_paths2 = []
    with open(read_path_file, 'r') as in1:
        for line in in1:
            read_paths.append(line)

    if read_path_file2 != "-99": # Determines if the paired end 
        with open(read_path_file2, 'r') as in2:
            for line in in2:
                read_paths2.append(line)

    return read_paths, read_paths2


def build_command(readpath, outdir, cutadapt_path, sickle_path,readpath2 = "-99"):
    """Return a cutadapt command-like string given the sample file(s)"""

    cmd = [cutadapt_path] # Begin building command.  Command components are stored as separate elements in a list, which will later be joined in a space-delimited fashion.

    if readpath2 == "-99": 
        cmd.append("--interleaved") # Modify command if paired-data is interleaved
        readpath2 = ''

    cmd.append('-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -f fastq --quality-base=33 -o') # All files are trimmed with these same common options

    ofile = readpath.split("/")[-1].strip("\n") # Get just the base filename and strip off the new line character.
    ofile = ofile.replace('.fq.gz', '_R1_cutadapt.fq.gz') # replaces either file extension with same common cutadapt extension
    ofile = ofile.replace('.fastq.gz', '_R1_cutadapt.fq.gz')
    ofile2 = ofile.replace('R1', 'R2')

    cmd.append(outdir + ofile)
    cmd.append('-p')
    cmd.append(outdir + ofile2)

    # And then append the read path
    cmd.append(readpath.strip("\n").replace(" ", "\ "))
    cmd.append(readpath2.strip("\n").replace(" ", "\ "))

    sfile1 = ofile.replace(".fq.gz","_sickle.fq")
    sfile2 = ofile2.replace(".fq.gz","_sickle.fq")
    sfile3 = ofile.replace('_R1_cutadapt.fq.gz', '_singles_cutadapt_sickle.fq')
    # Append sickle command
    cmd.append(f"; {sickle_path} pe -f {outdir}{ofile} -r {outdir}{ofile2} -t sanger -o {outdir}{sfile1} -p {outdir}{sfile2} -s {outdir}{sfile3}")

    # Remove cutadapt output after sickle finishes
    cmd.append(f"&& (rm {outdir}/{ofile}; rm {outdir}/{ofile2}); gzip {outdir}/{sfile1}; gzip {outdir}/{sfile2}; gzip {outdir}/{sfile3}")
    # And return the command as a string
    return ' '.join(cmd)


if __name__ == "__main__":
    """Main function."""

    parser = argparse.ArgumentParser(description='This script generates commands to be used with for trimming adapters from paired-end read data using cutadapt')
    parser.add_argument('-r', type=str, metavar='read_path_file', help='REQUIRED:  Full path to file with a list of read paths, one per line.  IMPORTANT: if paired reads are in separate files (i.e. not interleaved), create two files: one for forward reads and one for reverse reads. Pass the file containing paths to reverse read files to -r2')
    parser.add_argument('-r2', type=str, metavar='read_path_file2', default="-99", help='file containing paths to reverse read files')
    parser.add_argument('-o', type=str, metavar='output_directory', default='/panfs/roc/scratch/pmonnaha/Maize/Reads/Cutadapt/', help='Directory where the output of cutadapt and sickle will be written')
    parser.add_argument('-c', type=str, metavar='path_to_cutadapt', default='~/.local/bin/cutadapt', help='path to latest version of cutadapt')
    parser.add_argument('-s', type=str, metavar='path_to_sickle', default='~/software/sickle/sickle')

    args = parser.parse_args()

    #check python version
    assert sys.version_info[0] != 2, print("This script requires python3.  Do 'module load python/3.6.3'")

    home = str(Path.home())
    cutadapt_path = args.c.replace("~", home)
    sickle_path = args.s.replace("~", home)

    assert shutil.which(cutadapt_path), print("Did not find the CutAdapt executable at", cutadapt_path, ".  Be sure to provide exact path to executable file.")
    assert shutil.which(sickle_path), print("Did not find the Sickle executable at", sickle_path, ".  Be sure to provide exact path to executable file.")

    if args.o.endswith("/") is False: # Add "/" to end of output path if it doesn't already exist so that path to output file is formatted correctly
        args.o += "/"

    r_paths, r2_paths = get_read_paths(args.r, args.r2) # r2 = reverse reads; r = forward reads
    r_paths.sort()

    if args.r2 == "-99": # If data is interleaved, set null value (-99) for all elements in the list containing reverse reads
        r2_paths = ["-99" for j in r_paths]

    r2_paths.sort()

    for i, fq in enumerate(r_paths): # Loop over files, build commands, and print to output.
        cmd = build_command(fq, args.o, cutadapt_path, sickle_path, r2_paths[i])
        print(cmd)

