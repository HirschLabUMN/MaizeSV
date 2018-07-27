#!/bin/bash
#PBS -l mem=3gb,nodes=1:ppn=1,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi
#PBS -o /path/to/stdout/directory/cutadapt.o
#PBS -e /path/to/stderr/directory/cutadapt.e

# Load modules
module load cutadapt

# MODIFY THESE PATHS BELOW AS WELL AS -o AND -e paths above
CMD_LIST="/path/to/file/containing/cutadapt_commands.txt"
SUCCESS_LIST="/path/to/file/to/store/successful/cutadapt_runs.txt"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST} | rev | cut -d " " -f 1 | rev >> "${SUCCESS_LIST}"