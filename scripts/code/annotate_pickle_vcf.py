import allel
import subprocess
import argparse
import os
import pdb


#requires bgzip

def padBedFile(bed_file, gene_buff):
    new_bed = bed_file.replace("bed","tmp.bed")
    new_file = open(new_bed, 'w')
    # Create a SEPARATE bed entry for <gene_buff> distance upstream and downstream of gene for annotating VCF.
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.split()
            start = int(line[1])
            end = int(line[2])
            if int(line[1]) - gene_buff < 0: new_start = 0
            else: new_start = int(line[1]) - gene_buff
            new_end = int(line[2]) + gene_buff


            upstream = f"{line[0]}\t{new_start}\t{start}\t{line[3]}-us-{int(gene_buff / 1000)}kb\n"
            gene = f"{line[0]}\t{start}\t{end}\t{line[3]}-gene-{int(gene_buff / 1000)}kb\n"
            downstream = f"{line[0]}\t{end}\t{new_end}\t{line[3]}-ds-{int(gene_buff / 1000)}kb\n"

            new_file.write(gene)
            new_file.write(upstream)
            new_file.write(downstream)

    return(new_bed)




def annotateVCF(merged_vcf, bed_annt_file, outfile, annt_dist, SURVIVOR_ant_path):

    cmd1 = f"{SURVIVOR_ant_path} -b {bed_annt_file} -i {merged_vcf} --anno_distance {annt_dist} -o {merged_vcf.replace('.vcf', '.tmp.vcf')}"
    cmd2 = f"""awk '/^##INFO=<ID=SVTYPE/ {{ printf("##INFO=<ID=overlapped_Annotations,Number=.,Type=String,Description=\\"Overlapped Annotations\\">\\n");}} {{print;}}' {merged_vcf.replace('.vcf', '.tmp.vcf')} | gzip > {outfile}; rm {merged_vcf.replace('.vcf', '.tmp.vcf')}"""

    pp1 = subprocess.Popen(cmd1, shell = True)
    out1, err1 = pp1.communicate()
    pp2 = subprocess.Popen(cmd2, shell = True)
    out2, err2 = pp2.communicate()
    return([[out1, err1], [out2, err2]])

def pickleVCF(annt_vcf, template_vcf = -9, vcftools_perl_folder = "/usr/local/bin/", chroms = [i for i in range(1, 11)]):
    #MUST USE fields='*' in order to parse info relevant for SVs.
    #MUST USE numbers={'variants/overlapped_Annotations': X} to store multiple annotations.  Find 

    shuff_vcf = annt_vcf.replace(".vcf.gz", ".ord.vcf.gz")
    if template_vcf == -9:
        cmd = f"export PERL5LIB={vcftools_perl_folder}; {vcftools_perl_folder}/vcf-sort {annt_vcf} | bgzip > {shuff_vcf}; tabix -p vcf {shuff_vcf}"
        
    else:
        cmd = f"export PERL5LIB={vcftools_perl_folder}; {vcftools_perl_folder}/vcf-shuffle-cols -t {template_vcf} {annt_vcf} | {vcftools_perl_folder}/vcf-sort | bgzip > {shuff_vcf}; tabix -p vcf {shuff_vcf}"

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out1, err1 = p.communicate()

    assert 'failed' not in str(err1), print("Tabix sorting failed for master VCF")
    #Convert vcf to numpy array

    for chrom in chroms: #create npz file per chromosome
        chrom_vcf = shuff_vcf.replace(".vcf.gz", f".{chrom}.vcf.gz")
        base_cmd = f"tabix -H {shuff_vcf} {chrom} > {chrom_vcf.replace('.vcf.gz', '.hd.vcf')}"
        # pdb.set_trace()
        for pre in ["", "chr", "chr0"]: # This handles the fact that chromosomes are named differently across references
            # pdb.set_trace()
            cmd = base_cmd + f"; tabix {shuff_vcf} {pre}{chrom} > {chrom_vcf.strip('.gz')}"
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out1, err1 = p.communicate()
            if os.stat(chrom_vcf.strip('.gz')).st_size != 0: break

        cmd = f"cat {chrom_vcf.replace('.vcf.gz', '.hd.vcf')} {chrom_vcf.strip('.gz')} | bgzip > {chrom_vcf}; rm {chrom_vcf.replace('.vcf.gz', '.hd.vcf')}; rm {chrom_vcf.strip('.gz')}"
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out1, err1 = p.communicate()

        #write npz file for current chromosome
        allel.vcf_to_npz(chrom_vcf, annt_vcf.replace(".vcf.gz", f"{chrom}"), fields='*', numbers={'variants/overlapped_Annotations': 5}, overwrite=True, types = {'calldata/GQ': 'i2'})
        os.remove(chrom_vcf)
    #Write out master npz file for all chroms just in case.  Not currently using this.
    allel.vcf_to_npz(shuff_vcf, annt_vcf.replace(".vcf.gz", ""), fields='*', numbers={'variants/overlapped_Annotations': 5}, overwrite=True, types = {'calldata/GQ': 'i2'})
    return()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='vcf_info_file', required=True, help='tab delimited file with each line containing (tab-delimited) all of the vcfs for a reference genome that you wish to be merged.  The first, second, and third columns are reserved for: 1.) a reference ID in common for all vcfs on the line, 2.) the full path to the reference fasta, and 3.) a bed file with gene locations.  Each following column should be a comma-delimited (no spaces) list that has the vcf path followed by an identifier for that vcf (e.g. GS.L')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-s', type=str, metavar='output_suffix', required=True, help='Output Suffix')
    parser.add_argument('-mp', type=str, metavar='mergeSVcallers_path', required=False, default='/Users/pmonnahan/Documents/Research/Maize/MaizeSV/software/mergeSVcallers/mergeSVcallers')
    parser.add_argument('-sp', type=str, metavar='SURVIVOR_ant_path', required=False, default='/Users/pmonnahan/Documents/Research/Maize/MaizeSV/software/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant')
    # parser.add_argument('-r', type=str, metavar='reference_fasta_path_file', required=True, help="file containing the reference genome fasta path for each set of VCFs that you are looking to merge.  Must be in same order as -f")
    parser.add_argument('-ad', type=str, metavar='annotation_distance', required=False, default='2000', help = "distance from SV to buffer for looking for gene overlap when annotating merged vcf with gene info")
    parser.add_argument('-md', type=str, metavar='Merge_bp_distance', required=False, default='100', help = "Merge SVs with both breakpoints N BP away")
    parser.add_argument('-mr', type=str, metavar='Merge_recip_overlap', required=False, default='0', help = "Reciprocal overlap also required for merging")
    parser.add_argument('-b', type=int, metavar='buffer', required=False, default=2000, help = "Adds this amount to either side of gene boundary for annotation of variants")
    # parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')
    args = parser.parse_args()

    # print("BE SURE THAT LUMPY RESULTS HAVE BEEN FILTERED TO ONLY KEEP THE 'RETAINED' VARIANTS")

    VCF_dict = {}
    with open(args.f, 'r') as vcf_file:
        for line in vcf_file:
            line = line.split()
            ref_id = line[0]
            VCF_dict[ref_id] = [line[1]] #holds the fasta path for this reference
            VCF_dict[ref_id].append(line[2]) #holds the bed annotation
            for entry in line[3:]:
                vcf, ID = entry.split(",")
                VCF_dict[ref_id].append([vcf, ID])

    for i, ref in enumerate(VCF_dict.keys()):
        vcf = VCF_dict[ref][2][0]
        anntfile = f"{args.o}/{ref}.annt.{args.s}.gz"
        bed_file = padBedFile(VCF_dict[ref][1], args.b)
        annotateVCF(vcf, bed_file, anntfile, args.ad, args.sp) 
        pickleVCF(anntfile, anntfile)
        os.remove(bed_file)
