#!/usr/bin/env python
"""Parses a series of FastQC reports and prints out command lines for adapter
removal with cutadapt. We do it this way because cutadapt runtime scales
with the number of adapters to remove, so we only want to trim the ones that
are actually present in the reads. Takes two arguments:
    1) Directory of unzipped FastQC reports
    2) File with a list of read paths, one per line
"""

import sys
import os
import pipes
import re
import subprocess
import argparse


def get_file_list(fastqc_dir):
    """Return full paths to the fastqc reports given a directory full of
    unzipped fastqc reports."""
    # generate the full absolute path to the supplied directory
    f_path = os.path.abspath(os.path.expanduser(fastqc_dir))
    # Build the file list. We find the sub directories unter the main fastqc
    # directory, and then give the path to the 'fastqc_data.txt' file under it
    fastqc_dirs = [
        os.path.join(f_path, d, 'fastqc_data.txt')
        for d
        in os.listdir(f_path)
        if os.path.isdir(os.path.join(f_path, d))]
    return fastqc_dirs


def get_contam(f_report):
    """Parse the FastQC report to identify the contaminating adapters. If none
    are present, just return the universal adapter. Else, return the universal
    adapter and any other contaminating adapters."""
    # Start off with the universal adapter. This will always be searched.
    contaminants = ['AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT']
    # iterate through the file. This wuld be a really great time for Perl's
    # paragraph mode... We identify the line that starts with
    # ">>Overrepresented sequences" and take data until we find
    # >>END_MODULE. All the stuff in between are the contaminants.
    with open(f_report, 'r') as f:
        in_overrep = False
        for line in f:
            if line.startswith('Filename'):
                fname = line.strip().split()[1]
            # If we hit >>Overrepresented, then we are in the right section
            if line.startswith('>>Overrepresented'):
                in_overrep = True
                continue
            # Once we hit >>END, then we are done with overrepresented seqs
            if line.startswith('>>END'):
                in_overrep = False
            # If we are in the overrepresented sequences section and have a line
            # that starts with the hash mark, it is just the header and we skip
            # it.
            if in_overrep and line.startswith('#'):
                continue
            elif in_overrep:
                cont = line.strip().split()[0]
                # Check if 'N' is in our contaminant. If it is, then we want to
                # skip this, because it causes problems with trimming.
                if 'N' in contaminants:
                    continue
                else:
                    contaminants.append(line.strip().split()[0])
            else:
                continue
    return fname, contaminants


def get_read_paths(read_paths):
    """Return a dictionary where the key is the basename of the read (because
    fastqc reports that) and the value is the quoted full path of the read
    on MSI"""
    reads = {}
    with open(read_paths, 'r') as f:
        for line in f:
            p = pipes.quote(line.strip())
            f = os.path.basename(line.strip())
            reads[f] = p
    return reads


def build_command(f_name, contaminants, read_paths, outdir, interleaved):
    """Return a cutadapt command-like string given the sample file and the
    list of contaminants that are present in it."""
    # Get the proper path to the file
    path = read_paths[f_name]
    # We want to extract the barcode out of the filename. We do this by using a
    # regular expression.
    # barcode = re.search('[ATCG]{6}', f_name).group()

    cmd = 'zcat ' + path + " | head -40 | grep ' ' - " # take first 40 lines (first 10 reads) from fastq and grab out the read ID line which contains the index sequence 
    pp = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
    output = pp.communicate()

     # isolate index sequence from read ID line...the final split call will split the two indexes for a dual indexed fastq 
    lines = str(output[0]).split('\n') 

    for i,line in enumerate(lines[:-1]):
        temp_barcodes = str(line).split()[1].split(":")[-1].split("+")
        if len(temp_barcodes) == 2:

            # if "N" not in temp_barcodes[0] and "N" not in temp_barcodes[1]:
            barcodes = ['GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' + temp_barcodes[0].strip("N") + 'ATCTCGTATGCCGTCTTCTGCTTG', 'AATGATACGGCGACCACCGAGATCTACAC' + temp_barcodes[1].strip("N") + 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT']
            if i == 0:
                og_barcodes = barcodes
            elif og_barcodes != barcodes:
                print str("Errors in barcodes for " + f_name + ". Double check that command has correct barcode")
                # break
        elif len(temp_barcodes) == 1:
            # if "N" not in temp_barcodes[0]:
            barcodes = ['GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' + temp_barcodes[0].strip("N") + 'ATCTCGTATGCCGTCTTCTGCTTG']
            if i == 0:  
                og_barcodes = barcodes
            elif og_barcodes != barcodes:
                print "Errors in barcodes for " + f_name + ". Double check that command has correct barcode"
                # break

    # And start building the command string. We will not discard reads.
    cmd = ['cutadapt']
    # If the barcode is in the dictionary of index adapters, then we want to
    # add that to the list of adapters to trim
    for barcode in barcodes:
        cmd.append('-b')
        cmd.append(barcode)
    # Add the contaminants to trim
    for con in contaminants:
        cmd.append('-b')
        cmd.append(con)
    # Specify fastq format
    cmd.append('-f')
    cmd.append('fastq')
    # Specify 33-offset quality encoding
    cmd.append('--quality-base=33')
    # Specify an output file
    if f_name.endswith("fastq.gz"):
        ofile = f_name.replace('.fastq.gz', '_cutadapt.fastq.gz')
        if interleaved == 'false':
            f_name2 = path.replace('1.fastq.gz', '2.fastq.gz')
            ofile2 = f_name.replace('1.fastq.gz', '2_cutadapt.fastq.gz')
    elif f_name.endswith("fq.gz"):
        ofile = f_name.replace('.fq.gz', '_cutadapt.fastq.gz')
        if interleaved == 'false':
            f_name2 = path.replace('1.fq.gz', '2.fq.gz')
            ofile2 = f_name.replace('1.fq.gz', '2_cutadapt.fastq.gz')
    cmd.append('-o')
    cmd.append(outdir + ofile)
    if interleaved == 'false':
        cmd.append('-p')
        cmd.append(outdir + ofile2)
    # And then append the read path
    cmd.append(path)
    cmd.append(f_name2)

    # And return the command as a string
    return ' '.join(cmd)


if __name__ == "__main__":
    """Main function."""
    parser = argparse.ArgumentParser(description='This script generates commands to be used with for trimming adapters using cutadapt')

    parser.add_argument('-f', type=str, metavar='fastqc_dir', help='REQUIRED: Directory of unzipped FastQC reports')
    parser.add_argument('-r', type=str, metavar='read_path_file', help='REQUIRED: Full path to file with a list of read paths, one per line.  IMPORTANT: if paired reads are in separate files, only include the first file name')
    parser.add_argument('-i', type=str, metavar='interleaved', help='REQUIRED: are paired reads in separate files or interleaved: true or false')
    parser.add_argument('-o', type=str, metavar='output_directory', default='/panfs/roc/scratch/pmonnaha/Maize/Reads/Cutadapt/', help='Directory where the output of cutadapt will be written')
    parser.add_argument('-O', type=str, metavar='command_file', default='/home/hirschc1/pmonnaha/JobScripts/accessory/cutadapt_commands.txt', help='Full path to the file to which cutadapt commands will be written')

    args = parser.parse_args()

    outfile = open(args.O, 'w')

    if args.o.endswith("/") is False:
        args.o += "/"

    reports = get_file_list(args.f)
    r_paths = get_read_paths(args.r)
    for r in reports:
        fname, c = get_contam(r)
        if args.i == "false" and (fname.endswith("2.fq.gz") or fname.endswith("2.fastq.gz")):
            pass
        else:
            cmd = build_command(fname, c, r_paths, args.o, args.i)
            print cmd
            outfile.write(cmd + "\n")

