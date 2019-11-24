

import pickle
import numpy as np
import argparse
import pdb
import line_profiler
from itertools import compress
import os

def getInfo(NPZ_file, fields = ['variants/overlapped_Annotations', 'variants/CHROM', 'variants/ID', 'variants/POS', 'variants/END', 'variants/SVTYPE', 'calldata/GT', 'calldata/GQ', 'samples', 'variants/NSAMP', 'variants/SR', 'variants/PE']):
    '''Retrieve vector containing data in each of the specified fields'''
    info = []
    for f in fields:
        try: #if an expected field is not found in the npz/vcf file, a key error will result
            info.append(NPZ_file[f])
        except KeyError:
            print("Did not find field: ", f)
            info.append(-9)
    return(info)


def findOverlap(ref, alt_ref, NPZ, results_file, samples_file, gq_thresh, tmp_dir, min_ind, chrom_num = 10):
    unmatched_file = open(f"{tmp_dir}/{ref}-{alt_ref}_unmatched.txt", 'w') #Temporary file to hold variants that have not been matched for this pair of references.  We hold these in a file and wait to print them out because this variant may be matched in another reference
    matched = [] # Holds matched variant IDs which we will use to filter previously unmatched variants during final print step.
    homs = pickle.load(open(NPZ[ref][1], 'rb')) # Load dictionaries containing homologs for current reference

    results_output = open(results_file, 'a')
    samples_output = open(samples_file, 'a')

    for chrom in range(1, chrom_num+1): # Loop over chromosomes. Default is for 10 maize chromosomes
        print("Comparing Reference:", ref, "to Alternative Reference:", alt_ref, "for chromosome:", chrom)
        ref_npz = NPZ[ref][0].strip(".npz") + f"{chrom}.npz" #These npz files hold the annotated data from the VCFs.
        alt_npz = NPZ[alt_ref][0].strip(".npz") + f"{chrom}.npz"
        with np.load(ref_npz, allow_pickle=True) as vcf, np.load(alt_npz, allow_pickle=True) as alt_vcf: #Open files 
            Gene_ids, Chrom, num_vars, Pos, End, Svtype, GT, GQ, samples, Num_alt, SR, PE = getInfo(vcf) #Retrieve necessary info from vcfs
            alt_Genes, alt_Chrom, alt_vars, alt_Pos, alt_End, alt_Svtype, alt_GT, alt_GQ, alt_samples, alt_Num_alt, alt_SR, alt_PE = getInfo(alt_vcf)
            
            #This will catch instances where there are no variants for a particular chromosome
            if isinstance(alt_Genes, int): continue
            elif isinstance(Gene_ids, int): continue

            alt_Genes = list(alt_Genes)

            assert list(alt_samples) == list(samples), "Sample order is not the same across VCFs" #Sample order must be the same across VCFs in order to faithfully calculate distance between genotype matrices.

            for idx in range(0, np.shape(num_vars)[0]): # Loop over variants
                
                pos = Pos[idx]
                end = End[idx]
                svtype = Svtype[idx]
                length = str(int(end) - int(pos))
                gt = GT[idx]
                mgq = np.mean(GQ[idx]) #Mean genotype quality for current reference
                ref_idx = [k for k in range(0, len(samples)) if ref==samples[k]][0] #index of the genotype matrix corresponding to reference genotype
                gt_ref_2ref = sum(GT[idx][ref_idx]) #Genotype of current reference genotype mapped to the current reference
                gq_ref_2ref = GQ[idx][ref_idx] #Genotype quality current reference genotype mapped to the current reference
                gq_pass = [False if i < gq_thresh else True for i in GQ[idx]] # bool array to track whether genotype passes quality threshhold 

                #Get support unit (SU) info from Lumpy.  Two types are discordant paired reads (PE) and split reads (SR).
                num_alt, sr, pe = [-9, -9, -9]
                if hasattr(PE, "__len__"):
                    pe = PE[idx]
                    sr = SR[idx]
                    num_alt = Num_alt[idx]
                    
                gt_raw = np.array(list(compress(gt, gq_pass)))
                   
                raw_num_alt = np.sum(gt_raw) / 2
                raw_num_genos = np.shape(gt_raw)[0]

                num_hets = sum([1 for x in gt if sum(x)==1]) #Count number of heterozygous samples

                if GQ[idx][ref_idx] < gq_thresh: gt_ref_2ref = -2 #Is the current ref genotype of low quality? genotype is with respect to current ref not alt
                
                gene_ids = list(filter(None, Gene_ids[idx]))
                var_genes = len(gene_ids) #Number of genes that overlap with this variant

                if var_genes == 0:
                    # print(ref, "none", idx, var_genes, chrom, pos, end, -9, -9, -9, -9, -9, -9, -9, -9, svtype, -9, length, -9, -9, np.shape(gt)[0], -9, -9, -9, mgq, -9, gt_ref_2ref, -9, -9, -9, pe, sr, num_hets, num_alt)
                    unmatched_line = f"{ref} none {chrom}_{idx} {var_genes} {chrom} {pos} {end} -9 -9 -9 -9 -9 -9 -9 -9 {svtype} -9 {length} -9 -9 {raw_num_genos} -9 -9 -9 {mgq} -9 {gt_ref_2ref} -9 -9 -9 {gq_ref_2ref} -9 -9 -9 {pe} {sr} {num_hets} {raw_num_alt} {raw_num_genos} {raw_num_alt}\n"
                    unmatched_file.write(unmatched_line)
                else:        
                    for gene in gene_ids: # Loop over all genes that overlap with this variant
                        # pdb.set_trace()

                        # Determine location of overlap in current reference variant
                        try:
                            gene, location, buff = gene.split("-")
                            if ref == "PHB47": gene = gene.strip("m.g")
                        except ValueError:
                            continue

                        # Get homologs for current gene. This will throw KeyError if there are no homologs for this gene                       
                        try:
                            alt_homs = homs[gene][alt_ref] 
                            num_homs = len(alt_homs)
                        except KeyError:
                            #print UNMATCHED output::ISSUE:what if this variant is matched for another alt_ref. print all these to single file? Keep list of matched samples? and then cat unmatched variants that are not matched 
                            num_homs = 0
                            unmatched_line = f"{ref} none {chrom}_{idx} {var_genes} {chrom} {pos} {end} {gene} {location} {num_homs} -9 -9 -9 -9 -9 {svtype} -9 {length} -9 -9 {raw_num_genos} -9 -9 -9 {mgq} -9 {gt_ref_2ref} -9 -9 -9 {gq_ref_2ref} -9 -9 -9 {pe} {sr} {num_hets} {raw_num_alt} {raw_num_genos} {raw_num_alt}\n"
                            unmatched_file.write(unmatched_line)
                            continue

                        for i in range(0, np.shape(alt_vars)[0]): #Loop over variants in alt_vcf
                            alt_genes = [x.split("-")[0] for x in list(filter(None, alt_Genes[i]))] #Vector containing gene_IDs associated with this variant.  removes upstream or downstream annotation
                            alt_locs = [x.split("-")[1] for x in list(filter(None, alt_Genes[i]))] #Stores upstream or downstream annotation
                            
                            if alt_ref != "PHB47":
                                matches = [x for x,y in enumerate(alt_genes) if y in alt_homs] #Stores index of the gene of the alt_variant if the geneID is contained in the list of homologs associated with current gene in current ref
                            else:
                                matches = [x for x,y in enumerate(alt_genes) if y.strip("m.g") in alt_homs]

                            for match in matches: #Loop over alt_genes indexes that had a match.  Goal is to print a separate entry for every gene pair/variant match, e.g. if a single variant corresponds to multiple overlapped geneIDs 
                                alt_gene = alt_genes[match] #This will grab the different matched genes
                                alt_loc = alt_locs[match]
                                alt_pos = alt_Pos[i] #Variant info stays the same over this loop
                                alt_end = alt_End[i]
                                alt_type = alt_Svtype[i]
                                alt_len = str(int(alt_end) - int(alt_pos))
                                try:    
                                    alt_gt = alt_GT[i]
    
                                    alt_ref_idx = [k for k in range(0, len(samples)) if alt_ref == samples[k]][0] #Index of ALT ref genotype in genotype array corresponding CURRENT ref vcf

                                    alt_alt_idx = [k for k in range(0, len(alt_samples)) if alt_ref == alt_samples[k]][0] #Index of ALT ref genotype in genotype array corresponding ALT ref vcf
                                    ref_alt_idx = [k for k in range(0, len(alt_samples)) if ref == alt_samples[k]][0] #Index of CURRENT ref genotype in genotype array corresponding ALT ref vcf

                                    gt_alt_2ref = sum(GT[idx][alt_ref_idx]) #ALT ref genotype in genotype array corresponding CURRENT ref vcf

                                    gt_alt_2alt = sum(alt_GT[i][alt_alt_idx]) #Genotype of ALT ref genotype in genotype array corresponding ALT ref vcf
                                    gt_ref_2alt = sum(alt_GT[i][ref_alt_idx]) #CURRENT ref genotype in genotype array corresponding ALT ref vcf


                                    gq_alt_2ref = GQ[idx][alt_ref_idx] #Genotype qualities corresponding to above genotypes

                                    gq_alt_2alt = alt_GQ[i][alt_alt_idx]
                                    gq_ref_2alt = alt_GQ[i][ref_alt_idx]

                                    alt_gq_pass = [False if i < gq_thresh else True for i in alt_GQ[i]] #Boolian array of genotypes that pass quality threshold for the ALT V F
                                    dual_pass = [True if j * k == 1 else False for j,k in zip(gq_pass, alt_gq_pass)] # bool array for genotypes that pass in both

                                    gt_new = np.array(list(compress(gt, dual_pass))) #retrieve genotypes that pass genotype quality thresholds in both current and alt ref
                                    alt_gt_new = np.array(list(compress(alt_gt, dual_pass)))

                                    if GQ[idx][alt_ref_idx] < gq_thresh: gt_alt_2ref = -2 #Is the alt ref genotype of low quality? genotype is with respect to current ref not alt. set gt to -2 which is value for missing data

                                    if gt_alt_2ref == 2: #alt_ref genotype does not match reference
                                        alt_gt_new = np.where(alt_gt_new==0, 1, 0) #Flip genotypes of alt_ref. 

                                    if gt_ref_2ref == 2: #CURRENT ref does not match CURRENT ref.  These should ultimately be filtered out, but are left in to estimate their incidence.
                                        # gt_new = np.where(gt_new==0, 1, 0)
                                        pass

                                    
                                    amgq = np.mean(alt_GQ[i]) #mean genotype quality for ALT vcf

                                    flt_num_alt = np.sum(gt_new) / 2

                                    num_genos = np.shape(gt_new)[0] #Number of genotypes that remain after quality filtration

                                    if num_genos > min_ind and gt_alt_2ref != -2: #Only calculate distance if a sufficient number of individuals pass quality filtration and the alt_ref is not missing
                                        filt_dist = np.linalg.norm(gt_new-alt_gt_new) # Sum the differences in genotype scores between CURRENT and ALT variant                                  
                                        filt_max_dist = np.linalg.norm(np.zeros((num_genos, 2))-np.full((num_genos,2), 1)) #This calculates the maximum possible distance between arrays of size equal to filter genotyped array
                                        prop_dist_filt = filt_dist / filt_max_dist #genotype distance as a proportion of the maximum possible distance for this number of genotypes

                                        geno_diff = gt-alt_gt
                                        raw_dist = np.linalg.norm(geno_diff) #Distance calculated on unfiltered genotypes
                                        max_dist = np.linalg.norm(np.zeros((100, 2))-np.full((100,2), 1))

                                        prop_dist = raw_dist / max_dist #NOTE: This can end up being > 1 due to missing genotypes [-1, -1] being subtracted from an alt_genotype [1, 1] such that distance is 4 versus 2 (when genotypes are [0, 0] and [1, 1])

                                        #SOMEWHERE IN HERE CREATE FUNCTIONALITY TO IDENTIFY SAMPLES WITH DISCORDANT GENOTYPES
                                        #filt_dist = np.linalg.norm(gt_new-alt_gt_new) Nonzero entries in here should tell you mismatches
                                        
                                        filtered = [False if j * k == 1 else True for j,k in zip(gq_pass, alt_gq_pass)] 
                                        filt_samps = list([k for k in samples * filtered if k])
                                        Filt_samps = list(compress(samples, filtered))
                                        #pass_samps = np.array(list(compress(samples, dual_pass)))
                                        
                                        mismatches = [False if all(geno_diff[k]) == 0 else True for k in range(0, len(geno_diff))]
                                        mismatch_samples = list(compress(samples, mismatches))

                                        samples_line = f"{ref} {chrom} {chrom}_{idx} {alt_ref} {chrom} {chrom}_{i} {','.join(Filt_samps)} {','.join(mismatch_samples)}"
                                        samples_output.write(samples_line)
                                        # pdb.set_trace()
                                        #gt_new = np.array(list(compress(gt, dual_pass))) , should be able to use something like this to get sample names /indexes once you have mismatched indexes.
                                        #ALSO ADD NUMBER OF MISMATCHED GENOTYPES?
                                    else:
                                        filt_dist, filt_max_dist, prop_dist, prop_dist_filt, raw_dist, max_dist = [-9 for x in range(0,6)]
                                        
                                except KeyError:
                                    filt_dist = -9
                                    max_dist = -9
                                    prop_dist = -9
                                    mgq = -9
                                    amgq = -9
                                    print("no genotype info for one of the VCFs")

                                #print MATCHED output
                                matched_line = f"{ref} {alt_ref} {chrom}_{idx} {var_genes} {chrom} {pos} {end} {gene} {location} {num_homs} {chrom}_{i} {alt_pos} {alt_end} {alt_gene} {alt_loc} {svtype} {alt_type} {length} {alt_len} {filt_dist} {num_genos} {prop_dist_filt} {raw_dist} {prop_dist} {mgq} {amgq} {gt_ref_2ref} {gt_alt_2ref} {gt_ref_2alt} {gt_alt_2alt} {gq_ref_2ref} {gq_alt_2ref} {gq_ref_2alt} {gq_alt_2alt} {pe} {sr} {num_hets} {raw_num_alt} {raw_num_genos} {flt_num_alt}\n"
                                results_output.write(matched_line)
                                matched.append(f"{ref}-{chrom}_{idx}") #Hold info for matched variants to be filtered out when recovering previously unmatched variants.
    unmatched_file.close()
    results_output.close()
    samples_output.close()

    return(matched)

def recoverNoMatch(matched, results_file, tmp_dir):
    '''This function goes through all of the temporary files holding unmatched variants from each pairwise comparison.  It then prints the unmatched variant if it is not contained in the final, total collection of matched variants'''
    files = [f for f in os.listdir(tmp_dir) if "unmatched" in f] #Grab file names
    printed = []

    results_output = open(results_file, 'a')

    for file in files: 
        # pdb.set_trace()
        with open(f"{tmp_dir}/{file}", 'r') as unmatched:
            for line in unmatched:
                line1 = line.split()
                ID = f"{line1[0]}-{line1[2]}"
                if ID not in matched and ID not in printed:
                    results_output.write(line)
                    printed.append(ID)
        os.remove(f"{tmp_dir}/{file}")
    results_output.close()
    return()

            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the primary script that is used to cross-validate SV variants across reference genomes.  Users should have run the necessary scripts to prepare the VCFs (output of annotate_pickle_vcfs.py) and gene-keys (output of pickle_homologue_dicts.py).  The first output file (suffix: .results.txt) summarizes relavent info for each match, whereas the second output file (suffix: .samples.txt) contains the list of samples that were filtered and/or did not match across references for each match, which could be used for filtering sample genotypes when creating consensus calls across references.")
    parser.add_argument('-f', type=str, metavar='npz_file', required=True, help='tab delimited file with 3 colums: 1.)key name used in homologue dictionary (e.g. B73, W22, PHB47, and PH207) 2.) the npz file of the converted vcf (created with annotate_pickle_vcf.py), and 3.) the compressed dictionary containing homolog info (created by pickle_homologue_dicts.py)')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help="Directory to write the two output files")
    parser.add_argument('-s', type=str, metavar='output_prefix', default="SVmap")
    parser.add_argument('-gq', type=int, metavar='genotype_quality_threshold', default=0, help="Filter all genotypes below this threshold. Directly impacts matching and subsequently genotype distance calculations; Default = 0")
    parser.add_argument('-mi', type=int, metavar='minimum_individuals', default=0, help="minimum number of inidividuals (with sufficient genotype quality in both references) needed to calculate a distance between genotype matrices; Default = 0")
    parser.add_argument('-d', type=str, metavar='tmp_directory', default="/Users/pmonnahan/Documents")

    args = parser.parse_args()

    NPZ = {}

    out1 = f"{args.o}/{args.s}.gq{args.gq}.mi{args.mi}.results.txt"
    out2 = f"{args.o}/{args.s}.gq{args.gq}.mi{args.mi}.samples.txt"

    #write header of outputs
    with open(out1, 'w') as Out1:
        Out1.write("ref alt_ref varID var_genes chrom pos end gene location num_homs altID alt_pos alt_end alt_gene alt_loc svtype alt_type length alt_len filt_dist num_genos prop_dist_filt raw_dist prop_dist mgq amgq gt_ref_2ref gt_alt_2ref gt_ref_2alt gt_alt_2alt gq_ref_2ref gq_alt_2ref gq_ref_2alt gq_alt_2alt pe sr num_hets raw_num_alt raw_num_geno flt_num_alt\n")
    with open(out2, 'w') as Out2:
        Out2.write("ref chrom varID alt_ref alt_chrom altID FilteredSamples MisMatchedSamples\n")


    Matched = [] #This list will hold all matched variants from each pairwise comparison
    with open(args.f, 'r') as npz_file: #Parse input in file.  
        for line in npz_file:
            ref, npz, homs = line.strip().split() #Each line should have the reference identifier, npz (i.e. compressed vcf file), and compressed homolog dictiionary
            NPZ[ref] = [npz, homs]

    assert len(NPZ.keys()) > 1, print("Check format of npz_file.  Did not find files for multiple references.")


    for ref in NPZ:
        # pdb.set_trace()
        for alt_ref in NPZ:
            if ref != alt_ref: #Avoid comparing VCF to self
                matched = findOverlap(ref, alt_ref, NPZ, out1, out2, args.gq, args.d, args.mi) #Call main function to cross-validate a pair of VCFs
                Matched += matched
    Matched = list( dict.fromkeys(Matched) ) #Removes duplicates for faster searching
    recoverNoMatch(Matched, out1, args.d)
