#!/bin/bash
#PBS -l mem=3gb,nodes=1:ppn=1,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi
#PBS -o /path/to/stdout/directory/sambamba_merge.o
#PBS -e /path/to/stderr/directory/sambamba_merge.e


# MODIFY THESE PATHS BELOW AS WELL AS -o AND -e paths above
CMD_LIST="/path/to/file/containing/sambamba_commands.txt"
SUCCESS_LIST="/path/to/file/to/store/successful/sambamba_commands_runs.txt"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST} | rev | cut -d " " -f 1 | rev >> "${SUCCESS_LIST}"