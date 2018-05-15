#!/bin/bash
#PBS -l pmem=2500mb,nodes=1:ppn=12,walltime=72:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -o /home/hirschc1/pmonnaha/OandE/ComputeGenomeMask.o
#PBS -e /home/hirschc1/pmonnaha/OandE/ComputeGenomeMask.e
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi


module load parallel
module load java

REF_PATH="/home/hirschc1/pmonnaha/JobScripts/accessory/reference_paths.txt"
OUT="/home/hirschc1/pmonnaha/misc-files/gstrip/"

REFS=( $(awk '{print $2}' $REF_PATH) )

SV_CMD="/panfs/roc/msisoft/java/jdk1.8.0_144/bin/java -Xmx2g -cp /home/hirschc1/pmonnaha/software/svtoolkit/lib/SVToolkit.jar:/home/hirschc1/pmonnaha/software/svtoolkit/lib/gatk/GenomeAnalysisTK.jar org.broadinstitute.sv.apps.ComputeGenomeMask"

export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH}

for REF in "${REFS[@]}";do	

	NAME=$(basename $REF | cut -d "_" -f 1);
	CONS=($(grep ">" "$REF" | awk '{split($0,a," "); split(a[1],b,">"); print b[2]}'));

	for i in "${CONS[@]}";do echo "$SV_CMD -R "$REF" -O "${OUT}""${NAME}"."${i}".mask.fa -readLength 100 -sequence "$i""; $SV_CMD -R "$REF" -O "${OUT}""${NAME}"."${i}".mask.fa -readLength 100 -sequence "$i";
		done;

done | parallel --jobs 12 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD}


for f in "${REFS[@]}";do 
	NAME=$(basename $REF | cut -d "_" -f 1);
	cat "${OUT}""${NAME}".*.mask.fa >> "${OUT}""${NAME}".mask.fa | parallel --jobs 12 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD};
done