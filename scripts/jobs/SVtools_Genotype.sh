#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=24,walltime=24:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -o /PATH/TO/STDOUT/FILE.o 
#PBS -e /PATH/TO/STDERR/FILE.e
#PBS -M USER@umn.edu
#PBS -q mesabi

# Load required modules
module load parallel

# source the root file for CNVnator
source /panfs/roc/msisoft/root/6.06.06/bin/thisroot.sh

# MODIFY TO THE COMMANDS FILE GENERATED VIA Generate_SVtools_Genotype_commands.py, AS WELL AS THE STDOUT AND STDERR FILES ABOVE

CMDS="/path/to/svtools_genotype_command_file.txt"

# THESE COMMANDS CAN ALSO BE SPLIT INTO A SMALLER SET OF FILES AND RUN AS A TASK ARRAY.  UNCOMMENT LINE BELOW TO DO SO
#CMDS="/path/to/file/containing/svtools_genotype_commands_${PBS_ARRAYID}"

# Uncomment the next line for running parallel across multiple nodes
#parallel -v -v -t --jobs 24 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${TASKS}
# Uncomment the next line for single-node parallel
parallel < ${CMDS}
