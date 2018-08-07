#!/bin/bash -x

# The following import/export operations must also be placed within your ./bashrc file.  This ensures that the path/module information is properly inherited as subjobs are launched by the Queue.jar application
#Import necessary modules
module load java/jdk1.8.0_144
module load samtools
module load htslib/1.6
module load R/3.3.3
module load libdrmaa/1.0.13

#export paths for proper inheritance within Genome STRiP.  
SV_DIR="/home/hirschc1/pmonnaha/software/svtoolkit"
export LD_LIBRARY_PATH=${SV_DIR}:${LD_LIBRARY_PATH}
export SV_DIR
export PATH=${SV_DIR}:${PATH}
export LD_LIBRARY_PATH=/panfs/roc/msisoft/libdrmaa/1.0.13/lib/:${LD_LIBRARY_PATH}

# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 9 ]; then
cat << EOF

This script implements Genome STRiP's SVGenotyper Pipeline.  All arguments are required and must 
follow correct order.  User should also specify queue and requested resources via command line.  See note below for use of default values.

Usage: qsub GSTRiP_CNVDiscovery.sh -F \"VCFile refName outDir memGb bamListFile outName\"

VCFile: vcf that you would like ot genotype the samples in bamListFile

refName: e.g. W22, B73, PH207, or PHB47.  Paths are hardcoded.  Change if necessary.

outDir: Output Directory

memGb: number of Gb of RAM to request.  This should be ~2Gb less than requested from the HPC (e.g. qsub ... -l mem=...)

bamListFile: file containing list of bams to include in the analysis.  Must end with .list



EOF
  exit 1
fi


#Define variables
VCF="$1"
REF="$2"
OUT="$3"
MEM="$4"
BAM_LIST="$5"
NAME="$6"


# check that reference fasta exists
if ! [ -e "/home/hirschc1/pmonnaha/misc-files/gstrip/${REF}_chr1-10.fasta" ]; then
  echo "Did not find reference file: $1" >&2
  exit 1
fi

OE="${OUT}/OandE"

mkdir -p $OE
cd $OE

#Execute CNVDiscoveryPipeline
classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
java -Xmx${MEM}g -cp ${classpath} \
     org.broadinstitute.gatk.queue.QCommandLine \
     -S ${SV_DIR}/qscript/SVGenotyper.q \
     -S ${SV_DIR}/qscript/SVQScript.q \ 
     -cp ${classpath} \
     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
     -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
     -R /home/hirschc1/pmonnaha/misc-files/gstrip/${REF}_chr1-10.fasta \
     -I ${BAM_LIST} \
     -O ${OUT}/${NAME}
     -md /home/hirschc1/pmonnaha/misc-files/gstrip/${REF}_MetaData_E2-0 \
     -md /home/hirschc1/pmonnaha/misc-files/gstrip/${REF}_MetaData_E2-1 \
     -md /home/hirschc1/pmonnaha/misc-files/gstrip/${REF}_MetaData_E2-2 \
     -runDirectory ${OUT} \
     -jobLogDir ${OUT}/logs \
     -jobRunner Drmaa \
     -gatkJobRunner Drmaa \
     -P depth.parityCorrectionThreshold:null \
     -retry 10 \
     -parallelRecords 100 \
     -memLimit ${MEM} \
     -jobNative '-l walltime=12:00:00 -A hirschc1' \
     -run