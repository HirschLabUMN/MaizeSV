#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script prepares the data contained in multiple SV VCFs (corresponding to the same samples called against different reference genomes) for subsequent analysis with SVMap.py (which links the SVs across reference genomes).  
Requirements: SURVIVOR_ant, vcftools, and bgzip/tabix are all required for this script.
 Takes the following arguments:
    -f :  REQUIRED: tab delimited file with each line containing 1.) the reference genotype ID (e.g. B73; used for naming), 2.) a bed file with gene locations ONLY and 3.) the vcf file.')
    -o :  REQUIRED: Full path to output directory in which the pickled VCFs will be written
    -s :  REQUIRED: Output suffix. 
    -b : Buffer distance. Adds this amount to either side of gene boundary for annotation of variants.  Actually, creates a separate bed entry for upstream and downstream regiions and the associated gene ID, so that we can get specific info on whether variant overlapped with gene body and/or US/DS regions.
    -sp : full path to installation of SURVIVOR_ant.  Default is /Users/pmonnahan/Documents/Research/Maize/MaizeSV/software/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant.  See README for installation instructions of SURVIVOR_ant.
    -vp : path to the directory containing the perl modules of vcftools.  Default is that it is installed in /usr/local/bin/.
    -ad: Flag passed to SURVIVOR_ant specifying the distance from SV to buffer for looking for gene overlap when annotating vcf with gene info; this is accomodated by -b, which provides more explicit info regarding the location where the overlap is found, so this can be left at 0 (default).  

"""



import allel
import subprocess
import argparse
import os
import pdb


#requires bgzip

def padBedFile(bed_file, gene_buff):
    '''This function creates a new temporary bed file that has the original gene boundaries plus additional entries for the up and downstream regions, which are determined by a buffer specified by the -b flag'''
    new_bed = bed_file.replace("bed","tmp.bed")
    new_file = open(new_bed, 'w')
    # Create a SEPARATE bed entry for <gene_buff> distance upstream and downstream of gene for annotating VCF.
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.split()
            start = int(line[1])
            end = int(line[2])
            if int(line[1]) - gene_buff < 0: new_start = 0 #avoid creating a negative coordinate at beginning of chromosome
            else: new_start = int(line[1]) - gene_buff
            new_end = int(line[2]) + gene_buff

            #Format string entries corresponding to upstream, gene body, and downstream regions.  
            upstream = f"{line[0]}\t{new_start}\t{start}\t{line[3]}-us-{int(gene_buff / 1000)}kb\n"
            gene = f"{line[0]}\t{start}\t{end}\t{line[3]}-gene-{int(gene_buff / 1000)}kb\n"
            downstream = f"{line[0]}\t{end}\t{new_end}\t{line[3]}-ds-{int(gene_buff / 1000)}kb\n"

            #Write entries to file
            new_file.write(gene)
            new_file.write(upstream)
            new_file.write(downstream)

    return(new_bed)


def annotateVCF(merged_vcf, bed_annt_file, outfile, annt_dist, SURVIVOR_ant_path):
    '''Use SURVIVOR_ant software to annotate SV variants with entries created in padBedFile'''
    cmd1 = f"{SURVIVOR_ant_path} -b {bed_annt_file} -i {merged_vcf} --anno_distance {annt_dist} -o {merged_vcf.replace('.vcf', '.tmp.vcf')}" #SURVIVOR_ant command produces temporary output whose header is correctly formatted for subsequent steps via cmd2
    cmd2 = f"""awk '/^##INFO=<ID=SVTYPE/ {{ printf("##INFO=<ID=overlapped_Annotations,Number=.,Type=String,Description=\\"Overlapped Annotations\\">\\n");}} {{print;}}' {merged_vcf.replace('.vcf', '.tmp.vcf')} | awk '/^##FORMAT=<ID=GT,Number=1/ {{ printf("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\\"Genotype quality\\">\\n");}} {{print;}}' | grep -v ID=LENGTH | gzip > {outfile}; rm {merged_vcf.replace('.vcf', '.tmp.vcf')}""" #Command to fix header
    # pdb.set_trace()
    pp1 = subprocess.Popen(cmd1, shell = True) #Run cmd1
    out1, err1 = pp1.communicate() #Wait for it to finish
    pp2 = subprocess.Popen(cmd2, shell = True) #Run cmd2
    out2, err2 = pp2.communicate()
    return([[out1, err1], [out2, err2]])

def pickleVCF(annt_vcf, template_vcf = -9, vcftools_perl_folder = "/usr/local/bin/", chroms = [i for i in range(1, 11)]):
    '''This function takes the annotated VCF and: 1.) reorders columns, 2.) Parses by chromosome 3.) and create a compressed/pickled NPZ file containing per chromosome VCF data.'''

    #Reorder columns in the vcf based on ordering in the provided template vcfs.  This is necessary for SVMap to calculate genotypic distance correctly.
    shuff_vcf = annt_vcf.replace(".vcf.gz", ".ord.vcf.gz")
    if template_vcf == -9: 
        cmd = f"export PERL5LIB={vcftools_perl_folder}; {vcftools_perl_folder}/vcf-sort {annt_vcf} | bgzip > {shuff_vcf}; tabix -p vcf {shuff_vcf}"        
    else:
        cmd = f"export PERL5LIB={vcftools_perl_folder}; {vcftools_perl_folder}/vcf-shuffle-cols -t {template_vcf} {annt_vcf} | {vcftools_perl_folder}/vcf-sort | bgzip > {shuff_vcf}; tabix -p vcf {shuff_vcf}"

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out1, err1 = p.communicate()

    #Do not proceed with VCF conversion if reordering of columns failed.
    assert 'failed' not in str(err1), print("Tabix sorting failed for master VCF")

    #Convert vcf to numpy array
    for chrom in chroms: #create npz file per chromosome
        chrom_vcf = shuff_vcf.replace(".vcf.gz", f".{chrom}.vcf.gz")
        base_cmd = f"tabix -H {shuff_vcf} {chrom} > {chrom_vcf.replace('.vcf.gz', '.hd.vcf')}"
        # pdb.set_trace()
        for pre in ["", "chr", "chr0"]: # This handles the fact that chromosomes are named differently across references
            # pdb.set_trace()
            cmd = base_cmd + f"; tabix {shuff_vcf} {pre}{chrom} > {chrom_vcf.strip('.gz')}" #This will produce a file of size 0 if tabix cannot find the chromosome in the vcf
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out1, err1 = p.communicate()
            if os.stat(chrom_vcf.strip('.gz')).st_size != 0: break #The loop continues until the chromosome specification in the above tabix command is successful in retreiving the variants for that chromosome.

        cmd = f"cat {chrom_vcf.replace('.vcf.gz', '.hd.vcf')} {chrom_vcf.strip('.gz')} | bgzip > {chrom_vcf}; rm {chrom_vcf.replace('.vcf.gz', '.hd.vcf')}; rm {chrom_vcf.strip('.gz')}" #Readds original header to the chromosome vcf and removes intermediate fiiles
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) #execute above command
        out1, err1 = p.communicate()

        #write npz file for current chromosome
            #MUST USE fields='*' in order to parse info relevant for SVs.
            #MUST USE numbers={'variants/overlapped_Annotations': X} to store multiple annotations. 
        allel.vcf_to_npz(chrom_vcf, annt_vcf.replace(".vcf.gz", f"{chrom}"), fields='*', numbers={'variants/overlapped_Annotations': 5}, overwrite=True, types = {'calldata/GQ': 'i2'})
        os.remove(chrom_vcf)

    #Write out master npz file for all chroms just in case.  Not currently using this.
    allel.vcf_to_npz(shuff_vcf, annt_vcf.replace(".vcf.gz", ""), fields='*', numbers={'variants/overlapped_Annotations': 5}, overwrite=True, types = {'calldata/GQ': 'i2'})
    return()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script prepares the data contained in multiple SV VCFs (corresponding to the same samples called against different reference genomes) for subsequent analysis with SVMap.py (which links the SVs across reference genomes).\nRequirements: SUVIVOR_ant, vcftools, and bgzip/tabix')
    parser.add_argument('-f', type=str, metavar='vcf_info_file', required=True, help='tab delimited file with each line containing 1.) the reference genotype ID (e.g. B73; used for naming), 2.) a bed file with gene locations ONLY and 3.) the vcf file.')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-s', type=str, metavar='output_suffix', required=True, help='Output Suffix')
    parser.add_argument('-sp', type=str, metavar='SURVIVOR_ant_path', required=False, default='/Users/pmonnahan/Documents/Research/Maize/MaizeSV/software/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant')
    parser.add_argument('-vp', type=str, metavar='vcftools_perl_folder', required=False, default='/usr/local/bin/')
    parser.add_argument('-ad', type=str, metavar='annotation_distance', required=False, default='0', help = "distance from SV to buffer for looking for gene overlap when annotating merged vcf with gene info; this is accomodated by -b, which provides more explicit info regarding the location where the overlap is found, so this can be left at 0 (default)")
    parser.add_argument('-b', type=int, metavar='buffer', required=False, default=2000, help = "Adds this amount to either side of gene boundary for annotation of variants")
    args = parser.parse_args()

    #Loop over VCF files, annotate variants with geneIDs, and pickle the VCFs for fast retrieval in SVMap.py
    VCF_dict = {}
    with open(args.f, 'r') as vcf_file:
        for i, line in enumerate(vcf_file):
            line = line.split()
            ref = line[0]
            bed_file = line[1]
            vcf = line[2]
            if vcf.endswith("gz"):
                cmd = f"gunzip {vcf}"
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out1, err1 = p.communicate()
                vcf = vcf.strip(".gz")

            anntfile = f"{args.o}/{ref}.annt.{args.s}.vcf.gz"
            padded_bed_file = padBedFile(bed_file, args.b)
            annotateVCF(vcf, padded_bed_file, anntfile, args.ad, args.sp) 
            if i==0: template_vcf = anntfile #This is necessary for reordering sample names in the VCFs so that all are in same order across references.
            pickleVCF(anntfile, template_vcf, args.vp)
            os.remove(padded_bed_file)
