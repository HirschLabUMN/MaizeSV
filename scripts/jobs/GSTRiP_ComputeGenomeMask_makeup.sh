#!/bin/bash
#PBS -l mem=8gb,nodes=1:ppn=1,walltime=48:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi
#PBS -o /home/hirschc1/pmonnaha/OandE/gstrip_cgm_makeup.o
#PBS -e /home/hirschc1/pmonnaha/OandE/gstrip_cgm_makeup.e


# MODIFY THESE PATHS BELOW AS WELL AS -o AND -e paths above
CMD_LIST="/home/hirschc1/pmonnaha/junk/makeup_gstrip_commands.txt"

SV_DIR="/home/hirschc1/pmonnaha/software/svtoolkit"

export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH}

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

echo "${CMD}"
eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST}