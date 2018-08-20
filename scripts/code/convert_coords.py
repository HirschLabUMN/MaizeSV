from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import argparse
import pdb
import line_profiler
import subprocess
import csv
import os
import random
import string
import collections

def makeCoordDict(bedfile):
    coord_dict = {}
    with open(bedfile, 'r') as bed:
        for line in bed:
            line = line.strip("\n").split("\t")
            chrom = line[0]
            start = line[1]
            stop = line[2]
            geneID = line[3]
            try:
                coord_dict[chrom].append([start,stop,geneID])
            except KeyError:
                coord_dict[chrom] = [[start,stop,geneID]]
    # pdb.set_trace()
    return(coord_dict)


def getSequences(coord_dict, source_ref, temp_file, offset):
    with open(temp_file, 'w') as multi_fasta:
        for i, chrom in enumerate(SeqIO.parse(source_ref, "fasta")):
            try:
                for coord in coord_dict[chrom.id]:
                    seq1 = str(chrom.seq)[int(coord[0]):int(coord[0]) + offset]
                    seq2 = str(chrom.seq)[int(coord[1]) - offset:int(coord[1])]
                    Nstart = str(len(seq1) - seq1.count("N"))
                    Nstop = str(len(seq2) - seq2.count("N"))
                    multi_fasta.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2]}S;{Nstart}\n{seq1}\n")
                    multi_fasta.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2]}E;{Nstop}\n{seq2}\n")
            except KeyError: pass
    multi_fasta.close()
    return()

def getNewCoords(multi_fasta, Toblastdb, offset, thresshold):
    blastn_cline = NcbiblastnCommandline(query=multi_fasta, db=Toblastdb, evalue=0.001,
                                        outfmt=5, max_target_seqs=2, max_hsps=2, out=multi_fasta.replace(".fasta",".xml"))
    blastn_cline()
    # pdb.set_trace()
    results = NCBIXML.parse(open(multi_fasta.replace(".fasta",".xml")))
    newCoords = {}
    found = []
    genes = []
    for result in results:
        gene = result.query.split(";")[1][:-1]
        genes.append(gene)
        loc = result.query.split(";")[1][-1]
        qlen = int(result.query.split(";")[-1])
        chrom = result.query[0]
        for alignment in result.alignments: #Not searching for best hit
            for hsps in alignment.hsps:
                # pdb.set_trace()
                if abs(hsps.align_length - qlen) / qlen < thresshold:
                    if loc == "S":
                        found.append(gene)
                        pos = hsps.sbjct_start - hsps.query_start
                    elif loc == "E":
                        pos = hsps.sbjct_end + (qlen - hsps.query_end)
                    else: print("Did not find gene bound of query gene")
                    try:
                        newCoords[gene].append(int(pos))
                    except KeyError:
                        newCoords[gene] = [int(chrom), int(pos)]

    dups = [item for item, count in collections.Counter(found).items() if count > 1]
    if dups:
        print(f"duplicate entries found for: {','.join(dups)}")

    unfound = set(found) - set(genes)
    if unfound:
        print(f"did not find good matches for: {','.join(unfound)}")
    # Sort the entries in the dictionary prior to returning, such that the order will be chr, start, stop
    for ID in newCoords:
        if len(newCoords[ID]) != 3: print(f"Did not find complete coordinates for {ID}")
        newCoords[ID].sort()

    return(newCoords)

def verifyNewCoords(newCoords, bed_file, threshold):
    verCoords = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip("\n").split("\t")
            # pdb.set_trace()
            start = int(line[1])
            stop = int(line[2])
            length = stop - start
            geneID = line[3]
            nchrom, nstart, nstop = newCoords[geneID]
            nlen = nstop - nstart
            if abs(nlen - length) / length > threshold:
                print(f"bad match for {geneID}: old {start}-{stop} new {nstart}-{nstop}")
            else:
                verCoords[geneID] = newCoords[geneID]
    return(verCoords)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument('-s', type=str, metavar='bed_file_of_genes_to_convert', required=True, help='comma separated file where first entry is the parent/merged gene and all subsequent entries are corresponding genes from split model in alternative reference')
    parser.add_argument('-b', type=str, metavar='blast_folder', required=True, help='Full path directory containing chromosome-parsed blast databases')
    parser.add_argument('-r', type=str, metavar='source_ref', required=False, default = "/Users/pmonnahan/Documents/Research/MaizeSV/references/B73_genesOnly.bed", help='Full path to GSTRiP_RedundancyAnnotator.sh')
    parser.add_argument('-v', action="store_true")
    parser.add_argument('-t', type=str, metavar="tmp_dir", default=os.getcwd())
    parser.add_argument('-o', type=int, metavar="offset", default=2000, help="This determines the length of the sequence from the start and end positions to be used for blast to get alt ref coordinates")
    parser.add_argument('-T1', type=int, metavar="size_similarity_threshold", default=0.3, help="difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the offset value")
    parser.add_argument('-T2', type=int, metavar="size_similarity_threshold", default=0.7, help="difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the offset value")
    # parser.add_argument('-o', type=str, metavar='output_VCF_name', required=True, help = "Full path to final merged vcf")
    # parser.add_argument('-k', action='store_true', help="Keep temporary files")
    args = parser.parse_args()

    coord_dict = makeCoordDict(args.s)
    print("Made Coordinate Dictionary")
    tmp_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    tmp_file = f"{args.t}/tmp{tmp_str}.fasta" 
    getSequences(coord_dict, args.r, tmp_file, args.o)
    print("Retrieved sequences of coordinates")
    newCoords = getNewCoords(tmp_file, args.b, args.o, args.T1)
    verCoords = verifyNewCoords(newCoords, args.s, args.T2)
    for gene in verCoords:
        print("\t".join(verCoords[gene]) + "\t" + gene)
    os.remove(args.t)
    # match = getNewCoords(chrom, start, stop, args.r, args.b)

        

    # Must 
