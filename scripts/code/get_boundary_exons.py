'''
    File name: get_boundary_exons.py
    Author: Patrick Monnahan
    Date created: 09/01/18
    Python Version: 3.6
    Project: Split Genes
    Upstream of: convert_coords_wExons.py
    Downstream of: JMM's longest transcript code
    Description: This program is used to create the input for convert_coors_wExons.py.  It takes a list of gff annotations that have been filtered to include only a single transcript per gene and outputs bed files containing the first and last exons.
'''

import argparse
import pdb

def splitGenes(gene_file, gene_prefix_list):
    split_lols = [[] for pre in gene_prefix_list]
    with open(gene_file, 'r') as splits:
        for line in splits:
            for i, prefix in enumerate(gene_prefix_list):
                if line.startswith(prefix):
                    split_lols[i].append(line.strip())
    return(split_lols)

def makeBeds(annotation_list, splits, out_suffix_list, min_exon_size, split_size):
    for k, path in enumerate(annotation_list):
        j, m = (0, 0)
        with open(path, 'r') as gff:
            out = open(out_suffix_list[k] + str(j), 'w')
            p, n = (0, 0) # p tracks number of exons considered while n track number of exons passing size requirement
            for i, line in enumerate(gff):
                line = line.strip().split("\t")
                if line[0][0] == "#": continue
                gene = line[8].split("_")[0].split("=")[1]
                if gene in splits[k] and line[2] == "exon": # Loop over exons in annotation
                    p += 1
                    start = int(line[3])
                    stop = int(line[4])
                    gene = line[8].split("_")[0].split("=")[1]
                    if stop - start >= min_exon_size: 
                        n += 1
                        if "oldGene" not in locals(): # Check if oldGene has been initialized for this annotation (i.e. is this the first gene)
                            oldGene = gene
                            Lstart, Ustart = (start, start)
                            Lstop, Ustop = (stop, stop)
                        elif oldGene == gene: # Still looking at exons from same gene
                            if start < Lstart: # is current exon upstream of current most upstream
                                Lstart = start
                                Lstop = stop
                            elif stop > Ustop: #is current exon downstream of current most downstream
                                Ustart = start
                                Ustop = stop
                        else: # We have moved onto a new gene
                            if n / p < 0.6: print(f"gene {gene} only has {str(n)} out of {str(p)} exons of sufficient size {str(n/p)}") 
                            m += 1
                            p, n = (0, 0)
                            out.write("\t".join([line[0], str(Lstart), str(Lstop), oldGene + ".L"]) + "\n")
                            out.write("\t".join([line[0], str(Ustart), str(Ustop), oldGene + ".U"]) + "\n")
                            oldGene = gene
                            Lstart, Ustart = (start, start) # Reinitialize values for new gene
                            Lstop, Ustop = (stop, stop)
                            if split_size != -9: # To split output files or not
                                if m % split_size == 0:
                                    out.close()
                                    j += 1
                                    out = open(outs[k] + str(j), 'w')

            # WILL YOU MISS SOME GENES ENTIRELY 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "This program is used to create the input for convert_coors_wExons.py.  It takes a list of gff annotations that have been filtered to include only a single transcript per gene and outputs bed files containing the first and last exons.")
    parser.add_argument('-a', type=str, metavar='annotation_list', required=True, help='Comma separated string (no spaces) containing paths to annotation files.  annotation files should exons corresponding to a just a single transcript.')
    parser.add_argument('-s', type=str, metavar='query_genes', required=True, help='File with geneID (one per line) for which you wish to retrieve first and last exon')
    parser.add_argument('-p', type=str, metavar='prefix_list', required=True, help='Comma separated string with gene prefixes in SAME order as annotation list. e.g. Zm00001 for B73')
    parser.add_argument('-m', type=int, metavar="minimum_exon_size", default=30, help="Only use exons of this minimum size for blasting")
    parser.add_argument('-o', type=str, metavar='out_suffix_list', required=True, help="comma separated list of output prefixes in SAME order as annotation list")
    parser.add_argument('-v', action="store_true")
    parser.add_argument('-split', type=int, default=-9, help="split files into size of X genes (half number of entries) ")
    args = parser.parse_args()

    gffs = args.a.split(",")
    pres = args.p.split(",")
    outs = args.o.split(",")

    gene_ListofLists = splitGenes(args.s, pres)
    makeBeds(gffs, gene_ListofLists, outs, args.m, args.split)


