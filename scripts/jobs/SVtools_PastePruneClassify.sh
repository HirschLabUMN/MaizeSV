#!/bin/bash
#PBS -A hirschc1
#PBS -q mesabi

set -ex

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 6 ]; then
cat << EOF
Usage: qsub ${SCRIPT} -F \"vcf_dir merged.vcf out_prefix dist eval_param te.bed\"
This script implements svtools lsort and lmerge to combine variants across multiple 
individual VCFs generated with speedseq sv. All arguments are required and must 
follow correct order.

vcf_dir: input directory containing per-individual VCFs.  If you have called 
		SVs against multiple reference genomes, create separate 
		directories for each reference and call the script once per
		input directory

merged.vcf: the merged vcf that was used in svtools genotype/copynumber

out_prefix: prefix to use for output vcf.  Cannot already exist

dist:  max separation distance (bp) of adjacent loci in
           cluster [50]

eval_param: evaluating parameter for choosing best bedpe in a
            cluster(e.g. af=AlleleFrequency default:af)

te.bed: gzipped bed file from repeat masker or elsewhere that specifies
		te locations
EOF
  exit 1
fi

# check that input directory exists and is a directory
if ! [ -d "$1" ]; then
  echo "$1 is not a directory" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi

# Check that merged.vcf exists
if ! [ -e "$2" ]; then
  echo "$2 not found" >&2
  exit 1
fi

# Check that output files do not already exist
if [ -e "${3}.sv.gt.cn.vcf.gz" ]; then
  echo "${3}.sv.gt.cn.vcf.gz already exists" >&2
  exit 1
fi
if [ -e "${3}.sv.pruned.vcf.gz" ]; then
  echo "${3}.sv.pruned.vcf.gz already exists" >&2
  exit 1
fi
if ! [ -e "$6" ]; then
  echo "$6 does not exist" >&2
  exit 1
fi

VCF_DIR="$1"
VCF="$2"
OUT="$3"
DIST="$4"
EVAL="$5"
TE="$6"

# activate python 2.7 virtual environment
source activate py27

# Move into input directory and verify that gt and cn subdirectories exist
cd ${VCF_DIR}
if ! [ -e ./gt ]; then
	echo "gt directory not found" >&2
	exit 1
fi
if ! [ -e ./cn ]; then
	echo "cn directory not found" >&2
	exit 1
fi

# Create a file that lists the vcf files that resulted from svtools copynumber
ls -1 cn/*vcf > cn.list

# Combine genotyped vcfs together
svtools vcfpaste \
-m  ${VCF} \
-f cn.list \
-q \
| bgzip -c \
> ${OUT}.sv.gt.cn.vcf.gz

# If vcfpaste completes successfully, use svtools prune to filter out additional variants deemed to be identical
if [ $? -ne 0 ]; then
	zcat ${OUT}.sv.gt.cn.vcf.gz \
	| svtools afreq \
	| svtools vcftobedpe \
	| svtools bedpesort \
	| svtools prune -s -d ${DIST} -e \"${EVAL}\" \
	| svtools bedpetovcf \
	| bgzip -c > ${OUT}.sv.pruned.vcf.gz
else
	echo "Pruning pipeline failed"
fi

find ./gt/ -name "*.vcf" -printf "%f\n" | cut -d "_" -f 1 | awk '{print $1"\t2"}' > sex.txt

if [ $? -ne 0 ]; then
	zcat ${OUT}.sv.pruned.vcf.gz \
	| svtools classify \
 	-g sex.txt \
 	-a ${TE} \
 	-m large_sample \
	| bgzip -c > ${OUT}.ls.vcf.gz
fi

#Filter out BND calls
if [ $? -ne 0 ]; then
    zgrep -v BND ${OUT}.ls.vcf.gz \
    #| vawk --header '{if(I$RETAINED==1) print $0}' | bgzip -c > ${OUT}.ls.RT.vcf.gz
fi

