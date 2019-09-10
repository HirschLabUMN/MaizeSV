import pickle
import numpy as np
import argparse
import pdb
import line_profiler
from itertools import compress


# @profile
def findOverlap(ref, NPZ, gq_thresh = 50):
    homs = pickle.load(open(NPZ[ref][1], 'rb')) # Load dictionaries containing homologs for current reference
    with np.load(NPZ[ref][0], allow_pickle=True) as vcf:
        for idx in range(0, np.shape(vcf['variants/ID'])[0]): # Loop over variants
            gene_ids = vcf['variants/overlapped_Annotations'][idx]
            gene_ids = list(filter(None, gene_ids))
            chrom = vcf['variants/CHROM'][idx]
            chrom = chrom.replace("chr0", "")
            chrom = chrom.replace("chr", "")
            for gene in gene_ids: # Loop over all genes that overlap with this variant
                # pdb.set_trace()
                try: homs[gene].items()
                except KeyError: continue   
                for alt_ref, alt_genes in homs[gene].items(): # Get homologs for current gene. This will throw KeyError if there are no homologs for this gene
                    num_genes = len(alt_genes)
                    # pdb.set_trace()
                    try: NPZ[alt_ref][0] #NPZ[alt_ref][0] is filename for this reference
                    except KeyError: continue
                    with np.load(NPZ[alt_ref][0].strip(".npz") + f"{chrom}.npz", allow_pickle=True) as alt_vcf:
                        alt_Genes = list(alt_vcf['variants/overlapped_Annotations'])
                        for j, alt_gene in enumerate(alt_genes): 
                            for i in range(0, np.shape(alt_vcf['variants/ID'])[0]): #Loop over variants in alt_vcf
                                if alt_gene in list(filter(None, alt_Genes[i])): # Check if homolog for current gene is annotated for current variant in alt_vcf
                                    pos = vcf['variants/POS'][idx]
                                    end = vcf['variants/END'][idx]
                                    svtype = vcf['variants/SVTYPE'][idx]
                                    alt_pos = alt_vcf['variants/POS'][i]
                                    alt_end = alt_vcf['variants/END'][i]
                                    alt_type = alt_vcf['variants/SVTYPE'][i]
                                    # varID = vcf['variants/ID'][idx]
                                    # altID = alt_vcf['variants/ID'][i]
                                    alt_len = str(int(alt_end) - int(alt_pos))
                                    length = str(int(end) - int(pos))
                                    try:
                                        assert alt_vcf['samples'].all() == vcf['samples'].all(), "Sample order is not the same across VCFs"
                                        gt = vcf['calldata/GT'][idx]
                                        alt_gt = alt_vcf['calldata/GT'][i]
                                        # pdb.set_trace()
                                        alt_ref_idx = [k for k in range(0, len(vcf['samples'])) if alt_ref==vcf['samples'][k]][0]
                                        ref_idx = [k for k in range(0, len(vcf['samples'])) if ref==vcf['samples'][k]][0]
                                        alt_ref_gt = sum(vcf['calldata/GT'][idx][alt_ref_idx])
                                        ref_gt = sum(vcf['calldata/GT'][idx][ref_idx])
                                        gq_pass = [False if i < gq_thresh else True for i in vcf['calldata/GQ'][idx]]
                                        alt_gq_pass = [False if i < gq_thresh else True for i in alt_vcf['calldata/GQ'][i]]
                                        dual_pass = [True if j * k == 1 else False for j,k in zip(gq_pass, alt_gq_pass)]
                                        gt_new = np.array(list(compress(gt, dual_pass)))
                                        alt_gt_new = np.array(list(compress(alt_gt, dual_pass)))

                                        if vcf['calldata/GQ'][idx][alt_ref_idx] < gq_thresh: alt_ref_gt = -2
                                        if vcf['calldata/GQ'][idx][ref_idx] < gq_thresh: ref_gt = -2

                                        flipped = 0
                                        ref_flipped = 0
                                        if alt_ref_gt == 2:
                                            if ref_gt == 2:
                                                gt_new = np.where(gt_new==0, 1, 0)
                                                ref_flipped = 1
                                            elif ref_gt == 1: ref_flipped = -1
                                            elif ref_gt == -2: ref_flipped = -2
                                            #Flip genotypes of alt_ref. 
                                            alt_gt_new = np.where(alt_gt_new==0, 1, 0)
                                            flipped = 1 
                                        elif alt_ref_gt == -2:
                                            if ref_gt == 2:
                                                gt_new = np.where(gt_new==0, 1, 0)
                                                ref_flipped = 1
                                            elif ref_gt == 1: ref_flipped = -1
                                            elif ref_gt == -2: ref_flipped = -2
                                            flipped = -2
                                        elif alt_ref_gt == 1:
                                            if ref_gt == 2:
                                                gt_new = np.where(gt_new==0, 1, 0)
                                                ref_flipped = 1
                                            elif ref_gt == 1: ref_flipped = -1
                                            elif ref_gt == -2: ref_flipped = -2
                                            flipped = -1
                                        elif alt_ref_gt == 0:
                                            if ref_gt == 2:
                                                gt_new = np.where(gt_new==0, 1, 0)
                                                ref_flipped = 1
                                            elif ref_gt == 1: ref_flipped = -1
                                            elif ref_gt == -2: ref_flipped = -2
                                            #Alt ref gt is het.

                                        mgq = np.mean(vcf['calldata/GQ'][idx])
                                        amgq = np.mean(alt_vcf['calldata/GQ'][i])

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
                                    try:
                                        alt_tags = alt_vcf['variants/UTAGS'][i]
                                        tags = vcf['variants/UTAGS'][idx]
                                        print(ref, alt_ref, idx, chrom, pos, end, gene, num_genes, i, alt_pos, alt_end, alt_gene, svtype, alt_type, length, alt_len, tags, alt_tags, filt_dist, num_genos, prop_dist_filt, raw_dist, prop_dist, mgq, amgq, flipped, ref_flipped)
                                    except KeyError:
                                        print(ref, alt_ref, idx, chrom, pos, end, gene, num_genes, i, alt_pos, alt_end, alt_gene, svtype, alt_type, length, alt_len, filt_dist, num_genos, prop_dist_filt, raw_dist, prop_dist, mgq, amgq, flipped, ref_flipped)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='npz_file', required=True, help='tab delimited file with key name used in homologue dictionary (e.g. B73, W22, PHB47, and PH207) followed by the npz file of the converted vcf')
    # parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')
    args = parser.parse_args()

    NPZ = {}

    #print header of output
    print("ref alt_ref varID chrom pos end gene num_genes altID alt_pos alt_end alt_gene svtype alt_type length alt_len softwares1 softwares2 filt_dist num_genos prop_dist_filt raw_dist prop_dist mgq amgq alt_flipped ref_flipped")

    with open(args.f, 'r') as npz_file:
        for line in npz_file:
            ref, npz, homs = line.strip().split()
            NPZ[ref] = [npz, homs]
    for ref in NPZ:
        # pdb.set_trace()
        findOverlap(ref, NPZ)


        # len(bvcf['variants/overlapped_Annotations'][1])

# aa=os.listdir("/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/")
# aa=[x for x in aa if x.endswith("annt.vcf")]
# for a in aa:
#     allel.vcf_to_npz("/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/" + a, "/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/" + a.replace(".vcf",".npz"), fields="*",numbers={'variants/overlapped_Annotations': 5},overwrite=True)

#>>> bvcf.keys()
#['samples', 'calldata/CN', 'calldata/CNF', 'calldata/CNL', 'calldata/CNP', 'calldata/CNQ', 'calldata/FT', 'calldata/GL', 'calldata/GP', 'calldata/GQ', 'calldata/GSPC', 'calldata/GT', 'calldata/PL', 'variants/ALT', 'variants/CHROM', 'variants/CIEND', 'variants/CIPOS', 'variants/END', 'variants/FILTER_PASS', 'variants/GCFRACTION', 'variants/GCLENGTH', 'variants/GLALTFREQ', 'variants/GLALTSUM', 'variants/GLHETSUM', 'variants/GLINBREEDINGCOEFF', 'variants/GLREFFREQ', 'variants/GLREFSUM', 'variants/GSCALLRATE', 'variants/GSCLUSTERSEP', 'variants/GSCLUSTERSEPWEIGHTEDMEAN', 'variants/GSCLUSTERSEPWEIGHTEDMEDIAN', 'variants/GSCNALLELES', 'variants/GSCNCATEGORY', 'variants/GSCNDIST', 'variants/GSCNMAX', 'variants/GSCNMIN', 'variants/GSCNQUAL', 'variants/GSDEPTHINTERVAL', 'variants/GSDUPLICATEOVERLAP', 'variants/GSDUPLICATES', 'variants/GSDUPLICATESCORE', 'variants/GSELENGTH', 'variants/GSEXPMEAN', 'variants/GSGMMWEIGHTS', 'variants/GSM1', 'variants/GSM2', 'variants/GSNNONREF', 'variants/GSNONVARSCORE', 'variants/GSNVARIANT', 'variants/HOMLEN', 'variants/HOMSEQ', 'variants/ID', 'variants/IMPRECISE', 'variants/NOVEL', 'variants/ORIGINALID', 'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/SVLEN', 'variants/SVTYPE', 'variants/is_snp', 'variants/numalt', 'variants/overlapped_Annotations', 'variants/svlen']#