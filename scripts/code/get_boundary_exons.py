'''
    File name: ConvertGenicCoords.py
    Author: Patrick Monnahan
    Date created: 09/01/18
    Python Version: 3.6
    Project: Split Genes
    Upstream of: calcVarRatios.R
    Downstream of: JMM's longest transcript code
    Description: This script takes a series of split/merged candidates and produces a single pseudo-annotation file containing coordinates of each gene along with their corresponding split/merged genes.
'''

import argparse
import pdb
import csv
from pyfaidx import Fasta
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def removeDups(Split_gene_dict):  
    newDict = {}      
    for key in Split_gene_dict: # Removes duplicate entries from traditional name key dictionary
        for i, val in enumerate(Split_gene_dict[key]):
            if i == 0:
                newDict[key] = [val]
            elif val not in newDict[key]:
                newDict[key].append(val)
    return(newDict)

def BLAST(query, subject, debug=False):
    seq1 = SeqRecord(Seq(query), id="seq1")
    seq2 = SeqRecord(Seq(subject), id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")
    # Run BLAST and parse the output as XML
    if len(query) < 200:
        output = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5, max_target_seqs=1, task="blastn-short")()[0]
    else:
        output = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5, max_target_seqs=1)()[0]
    result = NCBIXML.read(StringIO(output))
    qstart, qstop, sstart, sstop = (-9,-9,-9,-9)
    for Alignment in result.alignments:
        for i, hsp in enumerate(Alignment.hsps):
            if i == 0: # Only going to look at the first (best) match
                qstart = min(hsp.query_start, hsp.query_end) # Account for fact that alignment could have various orientations
                qstop = max(hsp.query_start, hsp.query_end)
                sstart = min(hsp.sbjct_start, hsp.sbjct_end)
                sstop = max(hsp.sbjct_start, hsp.sbjct_end)
                if debug:
                    pdb.set_trace()
            else: break
    return(qstart, qstop, sstart, sstop)

def makeSplitGeneKeyDict(gene_key_file, source_prefix, alt_prefix):
    with open(gene_key_file, 'r') as split_gene_key:
        Alt_Parents = defaultdict(list) #Parent is from the reference that you are converting to...Alignment must be done for these
        Src_Parents = defaultdict(list) #Parent is from reference that you are converting from...Just use upper and lower bounds of children from gff.  no alignment necessary
        for i, line1 in enumerate(csv.reader(split_gene_key, delimiter = ",")):
            parent = line1[0]
            children = sorted(list(filter(None, line1[1:]))) # Remove empty entries in the child entries and sort entries (helps for identifying duplicate entries)
            if parent.startswith(source_prefix) and children[0].startswith(alt_prefix): 
                Src_Parents[parent].append([parent] + children)
                for child in children:
                    Src_Parents[child].append([parent] + children)
            elif parent.startswith(alt_prefix) and children[0].startswith(source_prefix): 
                Alt_Parents[parent].append([parent] + children)
                for child in children:
                    Alt_Parents[child].append([parent] + children)
    Alt_Parents = removeDups(Alt_Parents)
    Src_Parents = removeDups(Src_Parents)
    return(Alt_Parents, Src_Parents)

def makeCoordDict(annotation):
    coord_dict = {} # This dictionary uses genes as keys and simply stores coordinates
    with open(annotation, 'r') as gff:
        for line in gff:
            line = line.strip("\n").split("\t")
            if line[0][0] == "#": continue
            elif line[2] == "gene":
                gene = line[8].split("_")[0].split("=")[1].split(";")[0].split(".")[0]
                chrom = line[0]
                start = int(line[3])
                stop = int(line[4])
                coord_dict[gene] = [chrom, start, stop]
    return(coord_dict)

def convertCoords(SP_Dict, AP_Dict, srcCoords, alt_annotation, srcRef, altRef, min_exon_size, threshold, buff=0):
    AP_ParCoords = defaultdict(dict)
    newCoords = defaultdict(dict)
    srcSeq = Fasta(srcRef)
    altSeq = Fasta(altRef)
    p, n = (0, 0)
    oldGene = -9
    with open(alt_annotation, 'r') as gff: # loop over annotation file that you are converting coordinates
        for i, line in enumerate(gff):
            line = line.strip().split("\t")
            if line[0][0] == "#": continue
            gene = line[8].split("_")[0].split("=")[1].split(".")[0]
            start = int(line[3])
            stop = int(line[4])
            if gene in AP_Dict and line[2] == "gene": #             
                for fam in AP_Dict[gene]: # Loop over all split/merge gene sets associated with this gene
                    parent = fam[0]
                    key = ",".join(fam)
                    for child in fam[1:]:
                        childInfo = srcCoords[child]                       
                        try: 
                            AP_ParCoords[parent][key] += childInfo
                        except KeyError:
                            AP_ParCoords[parent][key] = childInfo
                        newCoords[child][key] = childInfo
    with open(src_annotation, 'r') as gff:
        for i, line in enumerate(gff):
            line = line.strip().split("\t")
            if line[0][0] == "#": continue
            gene = line[8].split("_")[0].split("=")[1].split(".")[0]
            start = int(line[3])
            stop = int(line[4])
            if gene in SP_Dict and line[2] == "exon":
                p += 1
                if stop - start >= min_exon_size: 
                    n += 1
                    if oldGene == -9: # Check if oldGene has been initialized for this annotation (i.e. is this the first gene)
                        oldGene = gene
                        exon_coords = [[start, stop]]
                    elif oldGene == gene: # Still looking at exons from same gene
                        exon_coords.append([start, stop])
                    else: # We have moved onto a new gene
                        if n / p < 0.6: print(f"gene {gene} only has {str(n)} out of {str(p)} exons of sufficient size {str(n/p)}") 
                        for par in SP_Dict[oldGene]: # Loop over all split/merge gene sets associated with this gene
                            key = ",".join(par)
                            parInfo = srcCoords[par[0]]
                            pchrom = parInfo[0].lstrip("chr").lstrip("0")
                            pstart = parInfo[1] - 1 - buff
                            pstop =  parInfo[2] - 1 + buff
                            subject = str(srcSeq[pchrom][pstart : pstop])
                            # pdb.set_trace()
                            # try:
                            #     subject = str(srcSeq[pchrom][pstart : pstop])
                            # except ValueError:
                            #     pdb.set_trace()
                            d = 0
                            new_coords = []
                            for i, exon in enumerate(exon_coords):
                                query = str(altSeq[pchrom][exon[0] : exon[1]])

                                qstart, qstop, sstart, sstop = BLAST(query, subject)

                                if qstart == -9:
                                    print(f"No blast results for exon {i+1} of {gene} with {par[0]}")
                                if (qstop - qstart) / len(query) < threshold: 
                                    d += 1
                                else:
                                    newStart = pstart + sstart 
                                    newStop = pstart + sstop
                                    new_coords += [newStart, newStop]
                            if (len(new_coords)/2) / len(exon_coords) < 0.6: print(f"gene {oldGene} only has {len(new_coords)/2} out of {len(exon_coords)} exons with sufficient alignment. SbjctLen {len(subject)}") 
                            else:
                                newStart = min(new_coords) 
                                newStop = max(new_coords)
                                newCoords[oldGene][key] = [pchrom, newStart, newStop]
                                newCoords[par[0]][key] = parInfo
                        exon_coords = [[start, stop]]
                        oldGene = gene
        for par in AP_ParCoords:
            for fam in AP_ParCoords[par]:
                coords = [x for x in AP_ParCoords[par][fam] if type(x) is int]
                chroms = set([x for x in AP_ParCoords[par][fam] if type(x) is str])
                assert len(chroms) == 1
                newCoords[par][fam] = [list(chroms)[0], min(coords), max(coords)]
    return(newCoords)

def verifyCoords(newCoords, altCoords, srcCoords, altPrefix, srcPrefix, threshold):
    """This function checks whether the new gene lengths determined via blasting are within a given threshold of the old distance determined from the annotation"""
    verCoords = {}
    for gene, fams in newCoords.items():
        try:
            oinfo = altCoords[gene]
        except KeyError:
            oinfo = srcCoords[gene]
        ochrom = oinfo[0].lstrip('chr0').strip("chr")
        olength = abs(oinfo[1] - oinfo[2])
        for fam, info in fams.items():
            nchrom = info[0]
            nstart = info[1]
            nstop = info[2]
            nlength = abs(nstop - nstart)
            pdb.set_trace()
            if ochrom != nchrom: print(f"Different chromosomes for new and old {gene}")
            elif abs(nlength - olength) / min(olength, nlength) > threshold: # Difference in new and old length expressed as a percent of the old length
                print(f"bad match for {gene} from {fam}: old {oinfo[2]}-{oinfo[1]} {olength} new {nstart}-{nstop} {abs(nlength)}")
            else:
                if nstop - nstart < 0: # Write entry such that start precedes stop
                    verCoords[gene] = [str(nchrom), str(nstop), str(nstart), gene, fam, str(oinfo[1]), str(oinfo[2])]
                else:
                    verCoords[gene] = [str(nchrom), str(nstart), str(nstop), gene, fam, str(oinfo[1]), str(oinfo[2])]
                print(f"good match for {gene} from {fam}: old {oinfo[2]}-{oinfo[1]} {olength} new {nstart}-{nstop} {abs(nlength)}")
    return(verCoords)

def writeBed(verCoords, outfile):
    with open(outfile, 'w') as out:
        for gene in verCoords:
            entry = [str(x) for x in verCoords[gene]]
            out.write("\t".join(entry) + "\n")
    return()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "This program is used to create the input for convert_coors_wExons.py.  It takes a list of gff annotations that have been filtered to include only a single transcript per gene and outputs bed files containing the first and last exons.")
    parser.add_argument('-s', type=str, metavar='split-gene-key', required=True, help='comma separated file where first entry is the parent/merged gene and all subsequent entries are corresponding genes from split model in alternative reference')
    parser.add_argument('-sa', type=str, metavar='annotation_list', required=True, help='Comma separated string (no spaces) containing paths to annotation files.  annotation files should exons corresponding to a just a single transcript.')
    parser.add_argument('-sr', type=str, metavar='source_ref', required=True, help='Reference fasta corresponding to reference converting FROM as in bed file')
    parser.add_argument('-aa', type=str, metavar='alt_annotation', required=True, help='Full path to alternative reference whose coordinates you are trying to convert to.  If using --global, then you must run makeblastdb on this fasta prior to running')
    parser.add_argument('-ar', type=str, metavar='alt_ref', required=True, help='Reference fasta corresponding to reference converting TO')
    parser.add_argument('-sp', type=str, metavar='source_prefix', required=True, help='gene prefix for source reference (e.g. Zm00008)')
    parser.add_argument('-ap', type=str, metavar='alt_prefix', required=True, help='gene prefix for alt reference (e.g. Zm00008)')
    parser.add_argument('-m', type=int, metavar="minimum_exon_size", default=15, help="Only use exons of this minimum size for blasting")
    parser.add_argument('-B', type=int, metavar="buffer", default=2000, help="For local blast to expected gene, how much buffer to add to both sides of the sbjct sequence")
    parser.add_argument('-o', type=str, metavar='out_file', required=True, help="Path to outfile")
    parser.add_argument('-t1', type=float, metavar='threshold', default=0.75, help="proportion of query that must align")
    parser.add_argument('-t2', type=float, metavar="size_similarity_threshold", default=3.0, help="proportion difference between gene size as function of old and new coordinates")
    parser.add_argument('-v', action="store_true")
    parser.add_argument('-split', type=int, default=-9, help="split files into size of X genes (half number of entries) ")
    args = parser.parse_args()

    Alt_Parents, Src_Parents = makeSplitGeneKeyDict(args.s, args.sp, args.ap)
    srcCoords = makeCoordDict(args.sa)
    altCoords = makeCoordDict(args.aa)
    newCoords = convertCoords(Src_Parents, Alt_Parents, srcCoords, args.aa, args.sr, args.ar, args.m, args.t1, args.B)
    verCoords = verifyCoords(newCoords, altCoords, srcCoords, args.sp, args.ap, args.t2)

    writeBed(verCoords, args.o)
