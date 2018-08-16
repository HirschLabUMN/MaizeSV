#!/bin/bash
#PBS -l pmem=2500mb,nodes=1:ppn=8,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi
#PBS -o /home/hirschc1/pmonnaha/OandE/sambamba_merge.o
#PBS -e /home/hirschc1/pmonnaha/OandE/sambamba_merge.e


# MODIFY THESE PATHS BELOW AS WELL AS -o AND -e paths above
CMD_LIST="/home/hirschc1/pmonnaha/JobScripts/accessory/merge_commands_E2_R1.txt"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST}