module load bcftools
source activate py27 #Activate python 2.7 virtualenv


set -ex

SCRIPT=`basename $0`
if [ "$#" -ne 1 ]; then
cat << EOF
Usage: sh ${SCRIPT} VCF_file
This script formats a VCF file produced by Genome STRiP so that it can be genotyped using Lumpy/SVtools.
EOF
    exit 1
fi

if ! [ -e "$1" ]; then
        echo "$1 not found" >&2
        exit 1
fi

VCF="$1"

gunzip $VCF || true
~/software/vcftools-vcftools-cb8e254/src/perl/vcf-sort ${VCF%.vcf.gz}.vcf | bgzip > ${VCF%.vcf.gz}.sort.vcf.gz
tabix -p vcf ${VCF%.vcf.gz}.sort.vcf.gz
zcat ${VCF%.vcf.gz}.sort.vcf.gz \
    | sed 's/GSCNCATEGORY,Number=1,Type=Integer/GSCNCATEGORY,Number=1,Type=String/' \
    | bcftools view -G \
    | vawk --header '{for (i = 1; i <= 7; i++) printf("%s\t", $i); printf("SVTYPE=%s;POS=%s;SVLEN=%s;END=%s;CIPOS=-10,10;CIEND=-10,10;CIPOS95=0,0;CIEND95=0,0\n", I$GSCNCATEGORY,$2,I$GSELENGTH,I$END)}' \
    | sed '/CIPOS,Number/a ##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">' \
    | sed '/CIEND,Number/a ##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">' \
    | sed 's/MIXED/DEL/g' > ${VCF%.vcf.gz}.fmted.vcf

create_coordinates -i ${VCF%.vcf.gz}.fmted.vcf -o ${VCF%.vcf.gz}.fmted_coord.txt

gzip ${VCF%.vcf.gz}.fmted.vcf
gzip ${VCF%.vcf.gz}.vcf
rm ${VCF%.vcf.gz}.sort.vcf.gz
rm ${VCF%.vcf.gz}.sort.vcf.gz.tbi

exit