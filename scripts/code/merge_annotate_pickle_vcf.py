import allel
import subprocess
import argparse
import os
import pdb


#requires bgzip

def prepVCF(vcf, vcftools_perl_folder = "/usr/local/bin/"):

    if vcf.endswith(".gz"):
        out = vcf.replace('.vcf.gz','.sort.vcf.gz')
        cmd = f"export PERL5LIB={vcftools_perl_folder}; zcat {vcf} | {vcftools_perl_folder}/vcf-sort | bgzip > {out}; tabix -p vcf {out}"
    elif vcf.endswith("vcf"):
        out = vcf.replace('.vcf','.sort.vcf.gz')
        cmd = f"export PERL5LIB={vcftools_perl_folder}; cat {vcf} | {vcftools_perl_folder}/vcf-sort | bgzip > {out}; tabix -p vcf {out}"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out1, err1 = p.communicate()
    return(out)

def mergeVCFs(vcf_info, merged_file_name, merge_dist, merge_ovlp, mergeSVcallers_path):
    ref_path = vcf_info[0]
    vcf_string = ""
    id_string = ""
    for vcf in vcf_info[2:]:
        Pvcf = prepVCF(vcf[0])
        vcf_string += Pvcf + ","
        id_string += vcf[1] + ","

    cmd = f"{mergeSVcallers_path} -a {ref_path} -f {vcf_string.strip(',')} -t {id_string.strip(',')} -s {merge_dist} -r {merge_ovlp} > {merged_file_name}"
    pp = subprocess.Popen(cmd, shell = True)
    out, err = pp.communicate()

    return(out, err)

def annotateVCF(merged_vcf, bed_annt_file, outfile, annt_dist, SURVIVOR_ant_path):

    cmd1 = f"{SURVIVOR_ant_path} -b {bed_annt_file} -i {merged_vcf} --anno_distance {annt_dist} -o {merged_vcf.replace('.vcf', '.tmp.vcf')}"
    cmd2 = f"""awk '/^##INFO=<ID=SVTYPE/ {{ printf("##INFO=<ID=overlapped_Annotations,Number=.,Type=String,Description=\\"Overlapped Annotations\\">\\n");}} {{print;}}' {merged_vcf.replace('.vcf', '.tmp.vcf')} | gzip > {outfile}; rm {merged_vcf.replace('.vcf', '.tmp.vcf')}"""

    pp1 = subprocess.Popen(cmd1, shell = True)
    out1, err1 = pp1.communicate()
    pp2 = subprocess.Popen(cmd2, shell = True)
    out2, err2 = pp2.communicate()
    return([[out1, err1], [out2, err2]])

def pickleVCF(annt_vcf):
    #MUST USE fields='*' in order to parse info relevant for SVs.
    #MUST USE numbers={'variants/overlapped_Annotations': X} to store multiple annotations.  Find 

    #Convert vcf to numpy array
    allel.vcf_to_npz(annt_vcf, annt_vcf.replace(".vcf.gz", ""), fields='*', numbers={'variants/overlapped_Annotations': 5}, overwrite=True)
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

    for ref in VCF_dict.keys():
        mergefile = f"{args.o}/{ref}.merge.{args.s}"
        mergeVCFs(VCF_dict[ref], mergefile, args.md, args.mr, args.mp)
        anntfile = f"{args.o}/{ref}.annt.{args.s}.gz"
        annotateVCF(mergefile, VCF_dict[ref][1], anntfile, args.ad, args.sp) 
        pickleVCF(anntfile)
