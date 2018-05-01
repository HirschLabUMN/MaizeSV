#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=8,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -q batch
#PBS -o /path/to/stdout/directory/fastqc.o
#PBS -e /path/to/stderr/directory/fastqc.e

# Load modules
module load fastqc

# MODIFY THE PATH BELOW AS WELL AS -o AND -e paths above
CMD_LIST="/path/to/file/containing/fastqc_commands.txt"

# CMD_LIST can be generated with "find <fastqdir> -name "*fastq*" -o -name "*fq*" -print | xargs -I {} echo "fastqc --noextract -t <number_of_threads> -o <outdir>" {} > <command_file>"
CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD}