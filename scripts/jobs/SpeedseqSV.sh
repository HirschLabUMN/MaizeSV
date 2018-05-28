#!/bin/bash
#PBS -l pmem=2500mb,nodes=1:ppn=8,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi


# MODIFY THESE PATHS BELOW AS WELL AS -o AND -e paths above
CMD_LIST="/home/hirschc1/pmonnaha/JobScripts/accessory/speedseqSV_commands_E2_R1.txt"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

echo "${CMD}"
eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST}
