#!/bin/bash
#PBS -l mem=31gb,nodes=1:ppn=12,walltime=72:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -o /home/hirschc1/pmonnaha/OandE/RepeatMasker.o
#PBS -e /home/hirschc1/pmonnaha/OandE/RepeatMasker.e
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

# Load modules
module load repeatmasker/4.0.5

cd /home/hirschc1/pmonnaha/misc-files/

RepeatMasker -pa 12 -s -species maize ../references/B73_chr1-10.fasta 
RepeatMasker -pa 12 -s -species maize ../references/PHB47_chr1-10.fasta
RepeatMasker -pa 12 -s -species maize ../references/PH207_chr1-10.fasta
RepeatMasker -pa 12 -s -species maize ../references/W22_chr1-10.fasta 
