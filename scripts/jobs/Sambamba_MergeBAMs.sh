#!/bin/bash
#PBS -l mem=16gb,nodes=1:ppn=1,walltime=18:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi

# Format arguments
SCRIPT=`basename $0`
if [ "$#" -ne 1 ]; then
cat << EOF
Usage: qsub ${SCRIPT} -F \"<Command_File>\"
        The purpose of this script is to submit, as a task array, jobs for merging BAM files (as provided in Command_File).  The Command_File can be generated via the Generate_MergeBAMs_commands.py or Generate_MergeExistingBAMs_commands.py scripts in the scripts/code subdirectory.
        
EOF
exit 1
fi

if ! [ -e "$1" ]; then
    echo "Command file: ${1} not found" >&2
    exit 1
fi

CMD_LIST="$1"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST}