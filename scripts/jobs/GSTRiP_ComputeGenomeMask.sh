#!/bin/bash
#PBS -l pmem=2500mb,nodes=1:ppn=10,walltime=72:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -o /home/hirschc1/pmonnaha/OandE/GSTRiP_ComputeGenomeMask.o
#PBS -e /home/hirschc1/pmonnaha/OandE/GSTRiP_ComputeGenomeMask.e
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

# Load necessary modules
module load parallel

REF_PATH="/home/hirschc1/pmonnaha/JobScripts/accessory/reference_paths.txt" # reference paths...one per line. 
OUT="/home/hirschc1/pmonnaha/misc-files/gstrip/" # Metadata directory

REFS=( $(awk '{print $1}' $REF_PATH) ) # Get second column from reference_paths.txt which contain the path to reference files

for REF in "${REFS[@]}";do	# Loop over references
	# Set directory and commands for ComputeGenomeMask
	SV_DIR="/home/hirschc1/pmonnaha/software/svtoolkit"
	SV_CMD="/panfs/roc/msisoft/java/jdk1.8.0_144/bin/java -Xmx2g -cp /home/hirschc1/pmonnaha/software/svtoolkit/lib/SVToolkit.jar:/home/hirschc1/pmonnaha/software/svtoolkit/lib/gatk/GenomeAnalysisTK.jar org.broadinstitute.sv.apps.ComputeGenomeMask"
	
	#Pull out the name of the reference as well as the contig names
	NAME=$(basename $REF | cut -d "_" -f 1);
	CONS=($(grep ">" "$REF" | awk '{split($0,a," "); split(a[1],b,">"); print b[2]}'));

	# Loop over chromosomes of a reference and run ComputeGenomeMask on each.
	for i in "${CONS[@]}";do
		echo "export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH};$SV_CMD -R "$REF" -O "${OUT}""${NAME}"."${i}".mask.fa -readLength 150 -sequence "$i""; 
		done | parallel --jobs 10 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD}; # Pipe commands for the 10 chromosomes to parallel
done

# Combine all of the per-chromosome masks within a reference
for f in "${REFS[@]}";do 
	NAME=$(basename $REF | cut -d "_" -f 1);
	cat "${OUT}""${NAME}".*.mask.fa >> "${OUT}""${NAME}".mask.fa | parallel --jobs 10 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD};
done