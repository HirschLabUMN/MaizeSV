#!/bin/bash
#PBS -l mem=3gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi


# And the command list
SCRIPT=`basename $0`
if [ "$#" -ne 1 ]; then
cat << EOF
Usage: qsub ${SCRIPT} -F \"<Command_File>\"
        This script is used to submit CutAdapt jobs, as listed in Command_File, as a task array.
        
EOF
exit 1
fi

if ! [ -e "$1" ]; then
    echo "Command file:$1 not found" >&2
    exit 1
fi

CMD_LIST="$1"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST} | rev | cut -d " " -f 1 | rev >> "${SUCCESS_LIST}"