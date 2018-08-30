import os
import argparse
import pdb
import line_profiler
import subprocess
import multiprocessing
from multiprocessing.pool import ThreadPool
from pprint import pprint

# @profile
def launchRA(vcf1, vcf2, refname, outpath, memGB, dupOvlThrsh, RAshPath):
    cmd = ['sh', RAshPath, vcf1, refname, outpath, memGB, dupOvlThrsh, vcf2]
    print(" ".join(cmd))
    # pdb.set_trace()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return(out, err)

def multiprocess(vcf_list, ref, outpath, memGB, dupOvlThrsh, RAshPath):
    pool = ThreadPool(multiprocessing.cpu_count())
    results = []
    tmpfiles = []
    for i, vcf1 in enumerate(vcf_list):
        for j, vcf2 in enumerate(vcf_list):
            vcf1 = vcf1.strip("\n")
            vcf2 = vcf2.strip("\n")
            outfile = outpath.replace(".vcf", f".tmp{i}{j}.vcf")
            tmpfiles.append(outfile)
            cmd = (vcf1, vcf2, ref, outfile, memGB, dupOvlThrsh, RAshPath)
            results.append(pool.apply_async(launchRA, cmd))
    pool.close()
    pool.join()
    return(results, tmpfiles)

def concatVCFs(vcf_list, vcftools_perl_folder, outfile):
    vcf_str = ' '.join(vcf_list)
    cmd = f"export PERL5LIB={vcftools_perl_folder}; {vcftools_perl_folder}/vcf-concat -p {vcf_str} | {vcftools_perl_folder}/vcf-sort -c | uniq > {outfile}"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    return(out, err)

def vawk(invcf, outvcf):
    cmd = "source activate py27; vawk --header '{if(I$GSDUPLICATES==\"NA\") print$0}' " + invcf + " | uniq > " + outvcf
    # pdb.set_trace()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    return(out, err)

def chr10(concatvcf, refname, outpath, memGB, dupOvlThrsh, RAshPath, vcftools):

    cmd = "source activate py27; vawk --header '{if($1==\"chr10\") print$0}' " + concatvcf + " > " + outpath.replace("vcf", "chr10.vcf")
    # pdb.set_trace()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    out4, err4 = launchRA(outpath.replace("vcf","chr10.vcf"), outpath.replace("vcf","chr10.vcf"), refname, outpath.replace("vcf", "chr10b.vcf"), memGB, dupOvlThrsh, RAshPath)
    out2, err2 = vawk(outpath.replace("vcf", "chr10b.vcf"), outpath.replace("vcf", "chr10c.vcf"))
    out3, err3 = concatVCFs([outpath, outpath.replace("vcf", "chr10c.vcf")], vcftools, outpath.replace(".vcf","2.vcf"))
    return([[out,out2,out3,out4],[err, err2, err3, err4]])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "This code will iteratively launch Genome STRiP's RedundancyAnnotator to produce a single merged VCF for all files provided in the input list file.")
    parser.add_argument('-f', type=str, metavar='vcf_file_list', required=True, help='file with vcfs listed one per line that you wish to be merged')
    parser.add_argument('-p', type=str, metavar='path_to_sh', required=False, default = "/home/hirschc1/pmonnaha/JobScripts/GSTRiP_RedundancyAnnotator.sh", help='Full path to GSTRiP_RedundancyAnnotator.sh')
    parser.add_argument('-m', type=str, metavar='memGB', required=False, default = "4")
    parser.add_argument('-r', type=str, metavar='refName', required=True, help = "W22, PH207, PHB47, or B73")
    parser.add_argument('-d', type=str, metavar='duplicateOverlapThreshold', required=False, default = "0.5", help = "Only variants with overlap greater than this threshold will be annotated.")
    parser.add_argument('-o', type=str, metavar='output_VCF_name', required=True, help = "Full path to final merged vcf")
    parser.add_argument('-v', type=str, metavar='vcftools_perl_folder', required=False, default="/home/hirschc1/pmonnaha/software/vcftools-vcftools-cb8e254/src/perl/", help = "")
    parser.add_argument('-k', action='store_true', help="Keep temporary files")
    parser.add_argument('--chr10_separate', action='store_true', help="run RA on chromosome 10 separately.  necessary for ph207 because gatk doesnt like the sorting of vcf tools")
    args = parser.parse_args()

    if args.r in ["W22", "B73", "PHB47", "PH207"]:
        with open(args.f, 'r') as vcf_file:
            vcf_list = []
            for vcf in vcf_file:
                vcf_list.append(vcf)
    else:
        print("-r must be W22, B73, PHB47, or PH207")

    # Prepare temporary file names
    tmpvcf = args.o.replace(".vcf", ".concat.vcf")
    prevcf = args.o.replace(".vcf", ".pre.vcf")
    prevcf2 = args.o.replace(".vcf", ".pre2.vcf")
    prevcf3 = args.o.replace(".vcf", ".pre3.vcf")

    # Run
    results, tmpfiles = multiprocess(vcf_list, args.r, args.o, args.m, args.d, args.p)
    out1, err1 = concatVCFs(tmpfiles, args.v, tmpvcf) 

    # The next launch call is where PH207 is failing at chr10 (because GATK walker is expecting chr10 to precede chr09...why?)
    # Need to vawk --header
    out2, err2 = launchRA(tmpvcf, tmpvcf, args.r, prevcf, args.m, args.d, args.p)
    out3, err3 = vawk(prevcf, prevcf2)
    out4, err4 = launchRA(prevcf2, prevcf2, args.r, prevcf3, args.m, args.d, args.p) # This is probably unnecessary.  I noticed that RA isnt calling preferred sites for some comparisons and thought by running it again, it might catch the prior no-calls
    out5, err5 = vawk(prevcf3, args.o)

    if args.chr10_separate:
        oo, ee = chr10(tmpvcf, args.r, args.o, args.m, args.d, args.p, args.v)

    
    # print stdout and stderr of each process
    for result in results:
        out, err = result.get()
        print(f"out: {out} err: {err}")
    for result in [[out1, err1], [out2, err2], [out3, err3]]:
        print(f"out: {result[0]} err: {result[1]}")

    # Delete temporary files
    if not args.k:
        for file in tmpfiles:
            os.remove(file)
        os.remove(tmpvcf)
        os.remove(prevcf)
