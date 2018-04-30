#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=8,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -q batch
#PBS -o /home/hirschc1/pmonnaha/OandE/fastqc.o
#PBS -e /home/hirschc1/pmonnaha/OandE/fastqc.e

# Load modules
module load fastqc

# And the command list
CMD_LIST="/panfs/roc/groups/14/hirschc1/pmonnaha/JobScripts/accessory/fastqc_preset1_recovery.txt"

# taken from the MSI parallel page at
# https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel

parallel --jobs 24 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${CMD_LIST}
