import pickle
import numpy as np
import allel
import argparse
import pdb
import line_profiler


# @profile
def findOverlap(ref, NPZ):
    vcf = np.load(NPZ[ref][0])
    homs = pickle.load(open(NPZ[ref][1], 'rb'))
    for idx in range(0, np.shape(vcf['variants/ID'])[0]):
        gene_ids = vcf['variants/overlapped_Annotations'][idx]
        gene_ids = list(filter(None, gene_ids))
        for gene in gene_ids:
            # pdb.set_trace()
            try:
                for alt_ref, alt_genes in homs[gene].items():
                    alt_vcf = np.load(NPZ[alt_ref][0])
                    alt_Genes = list(alt_vcf['variants/overlapped_Annotations'])
                    for j, alt_gene in enumerate(alt_genes):
                        # pdb.set_trace()
                        for i in range(0, np.shape(alt_vcf['variants/ID'])[0]):
                            if alt_gene[0] in list(filter(None, alt_Genes[i])):
                                pos = vcf['variants/POS'][idx]
                                chrom = vcf['variants/CHROM'][idx]
                                chrom = chrom.replace("chr0", "")
                                chrom = chrom.replace("chr", "")
                                end = vcf['variants/END'][idx]
                                svtype = vcf['variants/SVTYPE'][idx]
                                varID = vcf['variants/ID'][idx]
                                alt_pos = alt_vcf['variants/POS'][i]
                                alt_chrom = alt_vcf['variants/CHROM'][i]
                                alt_chrom = alt_chrom.replace("chr0", "")
                                alt_chrom = alt_chrom.replace("chr", "")
                                alt_end = alt_vcf['variants/END'][i]
                                alt_type = alt_vcf['variants/SVTYPE'][i]
                                altID = alt_vcf['variants/ID'][i]
                                alt_len = str(int(alt_end) - int(alt_pos))
                                length = str(int(end) - int(pos))
                                try:
                                    alt_tags = list(filter(None, alt_vcf['variants/UTAGS'][i]))
                                    tags = list(filter(None, vcf['variants/UTAGS'][i]))
                                    print(ref, alt_ref, varID, chrom, pos, end, gene, altID, alt_chrom, alt_pos, alt_end, alt_gene[0], alt_gene[1], svtype, alt_type, length, alt_len, ",".join(tags), ",".join(alt_tags))
                                except KeyError:
                                    print(ref, alt_ref, varID, chrom, pos, end, gene, altID, alt_chrom, alt_pos, alt_end, alt_gene[0], alt_gene[1], svtype, alt_type, length, alt_len)
            except KeyError: pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='npz_file', required=True, help='tab delimited file with key name used in homologue dictionary (e.g. B73, W22, PHB47, and PH207) followed by the npz file of the converted vcf')
    # parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')
    args = parser.parse_args()

    NPZ = {}
    with open(args.f, 'r') as npz_file:
        for line in npz_file:
            ref, npz, homs = line.strip().split()
            NPZ[ref] = [npz, homs]
    for ref in NPZ:
        findOverlap(ref, NPZ)


        # len(bvcf['variants/overlapped_Annotations'][1])

# aa=os.listdir("/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/")
# aa=[x for x in aa if x.endswith("annt.vcf")]
# for a in aa:
#     allel.vcf_to_npz("/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/" + a, "/Users/pmonnahan/Documents/Research/MaizeSV/data/VCFs/" + a.replace(".vcf",".npz"), fields="*",numbers={'variants/overlapped_Annotations': 5},overwrite=True)

#>>> bvcf.keys()
#['samples', 'calldata/CN', 'calldata/CNF', 'calldata/CNL', 'calldata/CNP', 'calldata/CNQ', 'calldata/FT', 'calldata/GL', 'calldata/GP', 'calldata/GQ', 'calldata/GSPC', 'calldata/GT', 'calldata/PL', 'variants/ALT', 'variants/CHROM', 'variants/CIEND', 'variants/CIPOS', 'variants/END', 'variants/FILTER_PASS', 'variants/GCFRACTION', 'variants/GCLENGTH', 'variants/GLALTFREQ', 'variants/GLALTSUM', 'variants/GLHETSUM', 'variants/GLINBREEDINGCOEFF', 'variants/GLREFFREQ', 'variants/GLREFSUM', 'variants/GSCALLRATE', 'variants/GSCLUSTERSEP', 'variants/GSCLUSTERSEPWEIGHTEDMEAN', 'variants/GSCLUSTERSEPWEIGHTEDMEDIAN', 'variants/GSCNALLELES', 'variants/GSCNCATEGORY', 'variants/GSCNDIST', 'variants/GSCNMAX', 'variants/GSCNMIN', 'variants/GSCNQUAL', 'variants/GSDEPTHINTERVAL', 'variants/GSDUPLICATEOVERLAP', 'variants/GSDUPLICATES', 'variants/GSDUPLICATESCORE', 'variants/GSELENGTH', 'variants/GSEXPMEAN', 'variants/GSGMMWEIGHTS', 'variants/GSM1', 'variants/GSM2', 'variants/GSNNONREF', 'variants/GSNONVARSCORE', 'variants/GSNVARIANT', 'variants/HOMLEN', 'variants/HOMSEQ', 'variants/ID', 'variants/IMPRECISE', 'variants/NOVEL', 'variants/ORIGINALID', 'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/SVLEN', 'variants/SVTYPE', 'variants/is_snp', 'variants/numalt', 'variants/overlapped_Annotations', 'variants/svlen']#