#!/bin/bash
#PBS -l mem=10gb,nodes=1:ppn=1,walltime=06:00:00 
#PBS -m abe 
#PBS -M pmonnaha@umn.edu 
#PBS -q mesabi
#PBS -o /home/hirschc1/pmonnaha/OandE/Lumpy.o
#PBS -e /home/hirschc1/pmonnaha/OandE/Lumpy.e

# load modules
module load samtools/1.4
module load java/jdk1.8.0_45
module load liblzma
source activate py27

# set paths
INPUT_DIR="/panfs/roc/scratch/pmonnaha/Maize/widiv_bams/ReadIDsort"
BED_DIR="/home/hirschc1/pmonnaha/misc-files"
OUTPUT_DIR="/panfs/roc/scratch/pmonnaha/Maize/widiv_bams/Lumpy"
INP="/panfs/roc/scratch/pmonnaha/Maize/widiv_bams/Cleaned"


# Get a list of the input alignments
BAMS=($(find ${INPUT_DIR} -name '*IDsorted.bam' | sort))
DISCS=($(find ${INPUT_DIR} -name '*_disc.bam' | sort))  # Get BAMS containing discordant read pairs
SPLTS=($(find ${INPUT_DIR} -name '*_splitters.bam' | sort))  # Get BAMS containing split reads

# Use PBS_ARRAYID to get a sample to process
C_ALN=${BAMS[${PBS_ARRAYID}]}

# Start processing the sample
bname=$(basename ${C_ALN})
sample=$(echo ${bname} | cut -d '_' -f 1-2)
ref=$(echo ${bname} | cut -d '_' -f 2)

# NG in the output file extension refers to the fact that these analyses include non-genic regions
/home/hirschc1/pmonnaha/software/lumpy-sv/bin/lumpyexpress -k -v -P -B ${C_ALN} -S ${SPLTS[${PBS_ARRAYID}]} -D ${DISCS[${PBS_ARRAYID}]} -o ${OUTPUT_DIR}/${sample}.NG.vcf -T ${OUTPUT_DIR}/${sample}.NG.tmp

svtyper -i ${OUTPUT_DIR}/${sample}.NG.vcf -B ${INP}/${sample}_Merged.bam -l ${OUTPUT_DIR}/${sample}.NG.json > ${OUTPUT_DIR}/${sample}.NG.gt.vcf

/home/hirschc1/pmonnaha/software/lumpy-sv/bin/lumpyexpress -k -v -P -x ${BED_DIR}/${ref}.NonGenicMask.bed -B ${C_ALN} -S ${SPLTS[${PBS_ARRAYID}]} -D ${DISCS[${PBS_ARRAYID}]} -o ${OUTPUT_DIR}/${sample}.vcf -T ${OUTPUT_DIR}/${sample}.tmp

svtyper -i ${OUTPUT_DIR}/${sample}.vcf -B ${INP}/${sample}_Merged.bam -l ${OUTPUT_DIR}/${sample}.json > ${OUTPUT_DIR}/${sample}.gt.vcf