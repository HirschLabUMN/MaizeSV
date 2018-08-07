#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=6,walltime=96:00:00
#PBS -A hirschc1
#PBS -m ae
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

# Load modules
module load parallel
module load liblzma

source activate py27

if [ "$#" -ne 1 ]; then
cat << EOF
Usage: qsub ${SCRIPT} -F \"command_file\"
EOF
  exit 1
fi

TASKS="$1"

# Run them with a huge parallel command. this is taken from the MSI help page
# that has an example of how to do this:
# https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel

cat ${PBS_NODEFILE}

# Uncomment the next line for running parallel across multiple nodes
#parallel -v -v -t --jobs 24 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${TASKS}
# Uncomment the next line for single-node parallel
parallel < ${TASKS}