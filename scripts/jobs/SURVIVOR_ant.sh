#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=6,walltime=96:00:00
#PBS -A hirschc1
#PBS -m ae
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

source activate py27

if [ "$#" -ne 3 ]; then
cat << EOF
Usage: qsub ${SCRIPT} -F \"vcf_bed_list annot_distance survivor_path \"

vcf_bed_list: comma-delimited file containing the VCF to be annotated in the first column
and the bed file with which to annotate in the second column.  VCFs must be unzipped.

annot_distance: Will annotate the variant with the bed feature if the feature is within
				this distance to the variant.
EOF
  exit 1
fi


#VCF functionality has not been tested!!!
while IFS="," read -r VCF BED EXT; do
	if [[ -z "$EXT" ]]; then
		echo "${3} -b ${BED} -i ${VCF} --anno_distance ${2} -o ${VCF%.vcf}.${2}annt.vcf"
    	$3 -b $BED -i $VCF --anno_distance $2 -o ${VCF%.vcf}.${2}tmp.vcf
        awk '/^##INFO=<ID=SVTYPE/ { printf("##INFO=<ID=overlapped_Annotations,Number=.,Type=String,Description=\"Overlapped Annotations\">\n");} {print;}' ${VCF%.vcf}.${2}tmp.vcf > ${VCF%.vcf}.${2}annt.vcf 
        rm ${VCF%.vcf}.${2}tmp.vcf
    else
    	EXT1=($(echo $EXT | tr '\t' ','))
    	echo "${3} -b ${BED} -i ${VCF} --anno_distance ${2} -o ${VCF%.vcf}.${2}annt.vcf"
    	$3 -b $BED -i $VCF --anno_distance $2 -o ${VCF%.vcf}.${2}annt.vcf -v $EXT1
    fi

done < "$1"

