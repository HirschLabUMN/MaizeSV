import pickle
import numpy as np
import argparse
import pdb
import line_profiler
from itertools import compress

def getInfo(NPZ_file, fields = ['variants/overlapped_Annotations', 'variants/CHROM', 'variants/ID', 'variants/POS', 'variants/END', 'variants/SVTYPE', 'calldata/GT', 'calldata/GQ', 'samples', 'variants/NSAMP', 'variants/SR', 'variants/PE']):
    info = []
    for f in fields:
        try: 
            info.append(NPZ_file[f])
        except KeyError:
            print("Did not find field: ", f)
            info.append(-9)
    return(info)


# @profile
def findOverlap(ref, NPZ, gq_thresh = 50):
    homs = pickle.load(open(NPZ[ref][1], 'rb')) # Load dictionaries containing homologs for current reference
    with np.load(NPZ[ref][0], allow_pickle=True) as vcf:
        # pdb.set_trace()
        Gene_ids, Chrom, num_vars, Pos, End, Svtype, GT, GQ, samples, Num_alt, SR, PE = getInfo(vcf)
        for idx in range(0, np.shape(num_vars)[0]): # Loop over variants
            # pdb.set_trace()
            matched = False # track whether this variant has been matched across references via homolog key
            gene_ids = list(filter(None, Gene_ids[idx]))
            chrom = Chrom[idx].replace("chr0", "").replace("chr", "")
            pos = Pos[idx]
            end = End[idx]
            svtype = Svtype[idx]
            length = str(int(end) - int(pos))
            gt = GT[idx]
            mgq = np.mean(GQ[idx]) #Mean genotype quality for current reference
            ref_idx = [k for k in range(0, len(samples)) if ref==samples[k]][0]
            gt_ref_2ref = sum(GT[idx][ref_idx]) #Genotype of current reference genotype mapped to the current reference
            gq_pass = [False if i < gq_thresh else True for i in GQ[idx]] # bool array to track whether genotype passes quality threshhold 
            var_genes = len(gene_ids) #Number of genes that overlap with this variant

            #Get support unit (SU) info from Lumpy.  Two types are discordant paired reads (PE) and split reads (SR).
            num_alt, sr, pe = [-9, -9, -9]
            if hasattr(PE, "__len__"):
                pe = PE[idx]
                sr = SR[idx]
                num_alt = Num_alt[idx]
                

            num_hets = sum([1 for x in gt if sum(x)==1]) #Count number of heterozygous samples

            if GQ[idx][ref_idx] < gq_thresh: gt_ref_2ref = -2 #Is the current ref genotype of low quality? genotype is with respect to current ref not alt
                         
            for gene in gene_ids: # Loop over all genes that overlap with this variant
                # pdb.set_trace()
                try: homs[gene].items() #See if there are any other homologs for this gene across all of the references stored in homs.
                except KeyError: continue   
                for alt_ref, alt_genes in homs[gene].items(): # Get homologs for current gene. This will throw KeyError if there are no homologs for this gene
                    num_homs = len(alt_genes)
                    # pdb.set_trace()
                    try: NPZ[alt_ref][0] #NPZ[alt_ref][0] is filename for this reference
                    except KeyError: continue
                    with np.load(NPZ[alt_ref][0].strip(".npz") + f"{chrom}.npz", allow_pickle=True) as alt_vcf:
                        alt_Genes, alt_Chrom, alt_vars, alt_Pos, alt_End, alt_Svtype, alt_GT, alt_GQ, alt_samples, alt_Num_alt, alt_SR, alt_PE = getInfo(alt_vcf)
                        alt_Genes = list(alt_Genes)
                        for j, alt_gene in enumerate(alt_genes): 
                            for i in range(0, np.shape(alt_vars)[0]): #Loop over variants in alt_vcf
                                if alt_gene in list(filter(None, alt_Genes[i])): # Check if homolog for current gene is annotated for current variant in alt_vcf
                                    alt_pos = alt_Pos[i]
                                    alt_end = alt_End[i]
                                    alt_type = alt_Svtype[i]
                                    alt_len = str(int(alt_end) - int(alt_pos))
                                    try:
                                        assert list(alt_samples) == list(samples), "Sample order is not the same across VCFs"
                                        alt_gt = alt_GT[i]
        
                                        alt_ref_idx = [k for k in range(0, len(samples)) if alt_ref == samples[k]][0]
                                        ref_idx = [k for k in range(0, len(samples)) if ref == samples[k]][0]

                                        alt_alt_idx = [k for k in range(0, len(alt_samples)) if alt_ref == alt_samples[k]][0]
                                        ref_alt_idx = [k for k in range(0, len(alt_samples)) if ref == alt_samples[k]][0]

                                        gt_alt_2ref = sum(GT[idx][alt_ref_idx])

                                        gt_alt_2alt = sum(alt_GT[i][alt_alt_idx])
                                        gt_ref_2alt = sum(alt_GT[i][ref_alt_idx])

                                        alt_gq_pass = [False if i < gq_thresh else True for i in alt_GQ[i]]
                                        dual_pass = [True if j * k == 1 else False for j,k in zip(gq_pass, alt_gq_pass)] # bool array for genotypes that pass in both
                                        gt_new = np.array(list(compress(gt, dual_pass))) #retrieve genotypes that pass genotype quality thresholds in both current and alt ref
                                        alt_gt_new = np.array(list(compress(alt_gt, dual_pass)))

                                        if GQ[idx][alt_ref_idx] < gq_thresh: gt_alt_2ref = -2 #Is the alt ref genotype of low quality? genotype is with respect to current ref not alt. set gt to -2 which is value for missing data

                                        if gt_alt_2ref == 2:
                                            #Flip genotypes of alt_ref. 
                                            alt_gt_new = np.where(alt_gt_new==0, 1, 0)

                                        if gt_ref_2ref == 2:
                                            # gt_new = np.where(gt_new==0, 1, 0)
                                            pass

                                        
                                        amgq = np.mean(alt_GQ[i])

                                        filt_dist = np.linalg.norm(gt_new-alt_gt_new)
                                        num_genos = np.shape(gt_new)[0]
                                        filt_max_dist = np.linalg.norm(np.zeros((num_genos, 2))-np.full((num_genos,2), 1))
                                        prop_dist_filt = filt_dist / filt_max_dist

                                        raw_dist = np.linalg.norm(gt-alt_gt)
                                        max_dist = np.linalg.norm(np.zeros((100, 2))-np.full((100,2), 1))
                                        prop_dist = raw_dist / max_dist
                                        
                                    except KeyError:
                                        filt_dist = -9
                                        max_dist = -9
                                        prop_dist = -9
                                        mgq = -9
                                        amgq = -9
                                        print("no genotype info for one of the VCFs")
                                    # try: # WAS USED FOR mergeSVcallers
                                    #     alt_tags = alt_vcf['variants/UTAGS'][i]
                                    #     tags = vcf['variants/UTAGS'][idx]
                                    #     print(ref, alt_ref, idx, var_genes, chrom, pos, end, gene, num_homs, i, alt_pos, alt_end, alt_gene, svtype, alt_type, length, alt_len, tags, alt_tags, filt_dist, num_genos, prop_dist_filt, raw_dist, prop_dist, mgq, amgq, gt_ref_2ref, gt_alt_2ref, gt_ref_2alt, gt_alt_2alt, pe, sr, num_hets, num_alt)
                                    #     matched = True
                                    # except KeyError:
                                    matched = True
                                    print(ref, alt_ref, idx, var_genes, chrom, pos, end, gene, num_homs, i, alt_pos, alt_end, alt_gene, svtype, alt_type, length, alt_len, filt_dist, num_genos, prop_dist_filt, raw_dist, prop_dist, mgq, amgq, gt_ref_2ref, gt_alt_2ref, gt_ref_2alt, gt_alt_2alt, pe, sr, num_hets, num_alt)
            if matched is False:
                print(ref, "none", idx, var_genes, chrom, pos, end, -9, -9, -9, -9, -9, -9, svtype, -9, length, -9, -9, np.shape(gt)[0], -9, -9, -9, mgq, -9, gt_ref_2ref, -9, -9, -9, pe, sr, num_hets, num_alt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='npz_file', required=True, help='tab delimited file with key name used in homologue dictionary (e.g. B73, W22, PHB47, and PH207) followed by the npz file of the converted vcf')
    parser.add_argument('-gq', type=int, metavar='genotype_quality_threshold', default=30)
    # parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')
    args = parser.parse_args()

    NPZ = {}

    #print header of output
    print("ref alt_ref varID var_genes chrom pos end gene num_homs altID alt_pos alt_end alt_gene svtype alt_type length alt_len filt_dist num_genos prop_dist_filt raw_dist prop_dist mgq amgq gt_ref_2ref, gt_alt_2ref, gt_ref_2alt, gt_alt_2alt, pe, sr, num_hets, num_alt")

    with open(args.f, 'r') as npz_file:
        for line in npz_file:
            ref, npz, homs = line.strip().split()
            NPZ[ref] = [npz, homs]
    for ref in NPZ:
        # pdb.set_trace()
        findOverlap(ref, NPZ, args.gq)


        # len(bvcf['variants/overlapped_Annotations'][1])

# aa=os.listdir("/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/")
# aa=[x for x in aa if x.endswith("annt.vcf")]
# for a in aa:
#     allel.vcf_to_npz("/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/" + a, "/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/" + a.replace(".vcf",".npz"), fields="*",numbers={'variants/overlapped_Annotations': 5},overwrite=True)

#>>> bvcf.keys()
#['samples', 'calldata/CN', 'calldata/CNF', 'calldata/CNL', 'calldata/CNP', 'calldata/CNQ', 'calldata/FT', 'calldata/GL', 'calldata/GP', 'calldata/GQ', 'calldata/GSPC', 'calldata/GT', 'calldata/PL', 'variants/ALT', 'variants/CHROM', 'variants/CIEND', 'variants/CIPOS', 'variants/END', 'variants/FILTER_PASS', 'variants/GCFRACTION', 'variants/GCLENGTH', 'variants/GLALTFREQ', 'variants/GLALTSUM', 'variants/GLHETSUM', 'variants/GLINBREEDINGCOEFF', 'variants/GLREFFREQ', 'variants/GLREFSUM', 'variants/GSCALLRATE', 'variants/GSCLUSTERSEP', 'variants/GSCLUSTERSEPWEIGHTEDMEAN', 'variants/GSCLUSTERSEPWEIGHTEDMEDIAN', 'variants/GSCNALLELES', 'variants/GSCNCATEGORY', 'variants/GSCNDIST', 'variants/GSCNMAX', 'variants/GSCNMIN', 'variants/GSCNQUAL', 'variants/GSDEPTHINTERVAL', 'variants/GSDUPLICATEOVERLAP', 'variants/GSDUPLICATES', 'variants/GSDUPLICATESCORE', 'variants/GSELENGTH', 'variants/GSEXPMEAN', 'variants/GSGMMWEIGHTS', 'variants/GSM1', 'variants/GSM2', 'variants/GSNNONREF', 'variants/GSNONVARSCORE', 'variants/GSNVARIANT', 'variants/HOMLEN', 'variants/HOMSEQ', 'variants/ID', 'variants/IMPRECISE', 'variants/NOVEL', 'variants/ORIGINALID', 'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/SVLEN', 'variants/SVTYPE', 'variants/is_snp', 'variants/numalt', 'variants/overlapped_Annotations', 'variants/svlen']#