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
import csv
from skbio import DNA, alignment
from pyfaidx import Fasta
from collections import defaultdict
from Bio import pairwise2
import operator
from math import ceil
import copy
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

gap_open = -10

gap_extend = -0.5

def slidingWindow(sequence,winSize,step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""
 
    # Verify the inputs
    try: it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")
 
    # Pre-compute number of chunks to emit
    numOfChunks = ceil(((len(sequence)-winSize)/step)+1)
    # Do the work
    for i in range(0,numOfChunks*step,step):
        yield sequence[i:i+winSize]

def removeDups(Split_gene_dict):  
    newDict = {}      
    for key in Split_gene_dict: # Removes duplicate entries from traditional name key dictionary
        for i, val in enumerate(Split_gene_dict[key]):
            if i == 0:
                newDict[key] = [val]
            elif val not in newDict[key]:
                newDict[key].append(val)
    return(newDict)

def BLAST(query, subject):
    seq1 = SeqRecord(Seq(query),
                   id="seq1")
    seq2 = SeqRecord(Seq("subject"),
                       id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse the output as XML
    output = NcbiblastpCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:
            print '****Alignment****'
            print 'sequence:', alignment.title
            print 'length:', alignment.length
            print 'e value:', hsp.expect
            print hsp.query
            print hsp.match
            print hsp.sbjct

def parseBlastResults(blastResults, newCoords, thresshold, breakFirst=False):
    noHits =[]
    genes = []
    found = []
    for result in blastResults:
        noHits.append(result.query)
        gene = result.query.split(";")[1]
        genes.append(gene)
        qlen = int(result.query.split(";")[-1])
        chrom = result.query[0]
        for alignment in result.alignments: 
            for hsps in alignment.hsps:
                if abs(hsps.align_length - qlen) / qlen < thresshold:
                    found.append(result.query)
                    if gene[-1] == "S":
                        pos = hsps.sbjct_start - hsps.query_start + 1
                    elif gene[-1] == "E":
                        pos = hsps.sbjct_end + (qlen - hsps.query_end)
                    else: print("Did not find gene bound of query gene")
                    try:
                        newCoords[gene].append(int(pos))
                    except KeyError:
                        newCoords[gene] = [int(chrom), int(pos)]
                    if breakFirst: break
            if breakFirst: break

    noHits = list(set(noHits) - set(found))
    return(newCoords, noHits, found)

def addGene(Gene_Dict, parent, children):
    try: # Only works if parent gene is already present in dictionary
        Gene_Dict[parent].append([parent] + children)
    except KeyError:
        Gene_Dict[parent] = [[parent] + children]
    for child in children:
        try:
            Gene_Dict[child].append([parent] + children)
        except KeyError:
            Gene_Dict[child] = [[parent] + children]
    return(Gene_Dict)

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
    # pdb.set_trace()
    return(coord_dict)

def convertCoords(SP_Dict, AP_Dict, srcCoords, alt_annotation, srcRef, altRef, min_exon_size, threshold, buff=0):
    AP_ParCoords = defaultdict(dict)
    newCoords = defaultdict(dict)
    srcSeq = Fasta(srcRef)
    altSeq = Fasta(altRef)
    p, n = (0, 0)
    oldGene = -9
    with open(alt_annotation, 'r') as gff:
        for i, line in enumerate(gff):
            line = line.strip().split("\t")
            if line[0][0] == "#": continue
            gene = line[8].split("_")[0].split("=")[1].split(".")[0]
            start = int(line[3])
            stop = int(line[4])
            if gene in AP_Dict and line[2] == "gene": # Loop over exons in annotation                
                for fam in AP_Dict[gene]:
                    parent = fam[0]
                    key = ",".join(fam)
                    for child in fam[1:]:
                        childInfo = srcCoords[child]                       
                        try:
                            AP_ParCoords[parent][key] += childInfo
                        except KeyError:
                            AP_ParCoords[parent][key] = childInfo
                        newCoords[child][key] = childInfo
                    # pdb.set_trace()
                oldGene = -9
                # pdb.set_trace()
    with open(alt_annotation, 'r') as gff:
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
                    # if oldGene == 'Zm00008a000981': pdb.set_trace()
                    if oldGene == -9: # Check if oldGene has been initialized for this annotation (i.e. is this the first gene)
                        oldGene = gene
                        exon_coords = [[start, stop]]
                    elif oldGene == gene: # Still looking at exons from same gene
                        exon_coords.append([start, stop])
                    else: # We have moved onto a new gene
                        if n / p < 0.6: print(f"gene {gene} only has {str(n)} out of {str(p)} exons of sufficient size {str(n/p)}") 
                        for par in SP_Dict[oldGene]:
                            key = ",".join(par)
                            parInfo = srcCoords[par[0]]
                            pstart = parInfo[1] - 1 - buff
                            pstop =  parInfo[2] - 1 + buff
                            subject = DNA(str(srcSeq[parInfo[0]][pstart : pstop]))
                            d = 0
                            new_coords = []
                            for exon in exon_coords:
                                query = DNA(str(altSeq[parInfo[0]][exon[0] : exon[1]]))
                                exon_length = abs(exon[1] - exon[0])
                                # winSize = int(exon_length * 3)
                                # stepSize = int(ceil(winSize * 0.25))
                                winSize = 10000
                                stepSize = 5000
                                if len(subject) > winSize:
                                    chunks = slidingWindow(subject, winSize, stepSize)
                                else:
                                    chunks = slidingWindow(subject, len(subject), len(subject))
                                aln_len = 0
                                for i, chunk in enumerate(chunks):
                                    try:
                                        Alignment = alignment.local_pairwise_align_ssw(query, chunk)
                                        rc_Alignment = alignment.local_pairwise_align_ssw(query.complement(), chunk)
                                        length = abs(Alignment[2][0][1] - Alignment[2][0][0])
                                        rc_length = abs(rc_Alignment[2][0][1] - rc_Alignment[2][0][0])
                                        
                                        if length < rc_length:
                                            Alignment = copy.deepcopy(rc_Alignment)
                                            length = abs(Alignment[2][0][1] - Alignment[2][0][0])
                                        if length / len(query) > threshold:
                                            if aln_len < length:
                                                aln_len = length
                                                aln_qcoords = Alignment[2][0]
                                                aln_scoords = Alignment[2][1]
                                                offset = i * stepSize
                                                newStart = pstart + (aln_scoords[0] + offset - aln_qcoords[0]) + 1
                                                assert aln_scoords[1] > aln_scoords[0]
                                                dist_from_subject_end = len(subject) - (offset + aln_scoords[1])
                                                newStop = pstop - dist_from_subject_end + (len(query) - aln_qcoords[1])
                                        # if oldGene == 'Zm00008a000981': pdb.set_trace()    
                                    except IndexError:
                                        print(f"Alignment error for chunk {i+1} {gene}. query len {len(query)} chunk len {len(chunk)} subject len {len(subject)}")
                                    # if oldGene == 'Zm00008a000981':
                                    #     print(oldGene, i * stepSize, i * stepSize + winSize, length, rc_length, aln_len, len(query))
                                    # pdb.set_trace()
                                # pdb.set_trace()
                                if aln_len / len(query) < threshold: d += 1
                                else:
                                    new_coords += [newStart, newStop]
                            # pdb.set_trace()
                            if (len(new_coords)/2) / len(exon_coords) < 0.6: print(f"gene {oldGene} only has {len(new_coords)/2} out of {len(exon_coords)} exons with sufficient alignment. SbjctLen {len(subject)}") 
                            else:
                                # pdb.set_trace()
                                newStart = min(new_coords)
                                newStop = max(new_coords)
                                newCoords[oldGene][key] = [parInfo[0], newStart, newStop]
                                newCoords[par[0]][key] = parInfo
                        exon_coords = [[start, stop]]
                        oldGene = gene
        # pdb.set_trace()
        for par in AP_ParCoords:
            for fam in AP_ParCoords[par]:
                coords = [x for x in AP_ParCoords[par][fam] if type(x) is int]
                chroms = set([x for x in AP_ParCoords[par][fam] if type(x) is str])
                assert len(chroms) == 1
                # print("1,",newCoords[par][fam])
                newCoords[par][fam] = [list(chroms)[0], min(coords), max(coords)]
                # print("2,",newCoords[par][fam])
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
    # Coords = makeCoordDict(gffs)
    newCoords = convertCoords(Src_Parents, Alt_Parents, srcCoords, args.aa, args.sr, args.ar, args.m, args.t1, args.B)
    verCoords = verifyCoords(newCoords, altCoords, srcCoords, args.sp, args.ap, args.t2)
    pdb.set_trace()

    writeBed(verCoords, args.o)



    # def makeBeds(annotation_list, splits, parCoordDict, out_suffix_list, min_exon_size, split_size):
#     for k, path in enumerate(annotation_list):
#         j, m = (0, 0)
#         with open(path, 'r') as gff:
#             if split_size == -9: j = ".bed"
#             out = open(out_suffix_list[k] + str(j), 'w')
#             p, n = (0, 0) # p tracks number of exons considered while n track number of exons passing size requirement
#             for i, line in enumerate(gff):
#                 line = line.strip().split("\t")
#                 if line[0][0] == "#": continue
#                 gene = line[8].split("_")[0].split("=")[1]
#                 if gene in splits and line[2] == "exon": # Loop over exons in annotation
#                     p += 1
#                     start = int(line[3])
#                     stop = int(line[4])
#                     parent = splits[gene][0][0]
#                     if stop - start >= min_exon_size: 
#                         n += 1
#                         if "oldGene" not in locals(): # Check if oldGene has been initialized for this annotation (i.e. is this the first gene)
#                             oldGene = gene
#                             Lstart, Ustart = (start, start)
#                             Lstop, Ustop = (stop, stop)
#                         elif oldGene == gene: # Still looking at exons from same gene
#                             if start < Lstart: # is current exon upstream of current most upstream
#                                 Lstart = start
#                                 Lstop = stop
#                             elif stop > Ustop: #is current exon downstream of current most downstream
#                                 Ustart = start
#                                 Ustop = stop
#                         else: # We have moved onto a new gene.  Do alignments here
#                             if n / p < 0.6: print(f"gene {gene} only has {str(n)} out of {str(p)} exons of sufficient size {str(n/p)}") 
#                             m += 1
#                             p, n = (0, 0)
#                             out.write("\t".join([line[0], str(Lstart), str(Lstop), oldGene + ".L", parent, parCoordDict[parent][2], str(parCoordDict[parent][0]), str(parCoordDict[parent][1])]) + "\n")
#                             out.write("\t".join([line[0], str(Ustart), str(Ustop), oldGene + ".U", parent, parCoordDict[parent][2], str(parCoordDict[parent][0]), str(parCoordDict[parent][1])]) + "\n")
#                             oldGene = gene
#                             Lstart, Ustart = (start, start) # Reinitialize values for new gene
#                             Lstop, Ustop = (stop, stop)
#                             if split_size != -9: # To split output files or not
#                                 if m % split_size == 0:
#                                     out.close()
#                                     j += 1
#                                     out = open(outs[k] + str(j), 'w')
                        # par1:[[par1,geneA,geneB],[par1,geneC,geneD]]
                        # geneC: [[par1,geneC,geneD],[par2,geneC,geneE]]
            #             # Alt_Par_Coords[parent][prog]
            # if gene in SP_Dict and line[2] == "exon":
            #     p += 1
            #     if stop - start >= min_exon_size: 
            #         n += 1
            #         if "oldGene" not in locals(): # Check if oldGene has been initialized for this annotation (i.e. is this the first gene)
            #             oldGene = gene
            #             Lstart, Ustart = (start, start)
            #             Lstop, Ustop = (stop, stop)
            #         elif oldGene == gene: # Still looking at exons from same gene
            #             if start < Lstart: # is current exon upstream of current most upstream
            #                 Lstart = start
            #                 Lstop = stop
            #             elif stop > Ustop: #is current exon downstream of current most downstream
            #                 Ustart = start
            #                 Ustop = stop
            #         else: # We have moved onto a new gene.  Do alignments here
            #             if n / p < 0.6: print(f"gene {gene} only has {str(n)} out of {str(p)} exons of sufficient size {str(n/p)}") 
            #             p, n = (0, 0)
                        
            #             for par in SP_Dict[gene]:
            #                 key = ",".join(par)
            #                 parInfo = srcCoords[par[0]]
            #                 pstart = parInfo[1] - 1 - buff
            #                 pstop =  parInfo[2] - 1 + buff
            #                 subject = DNA(str(srcSeq[parInfo[0]][pstart : pstop]))
            #                 # get seqs from Lstart: Lstop and Ustart:Ustop and align
            #                 lower_query = DNA(str(altSeq[chrom][Lstart - 1 : Lstop - 1])) # Pull sequence from reference using current gene coordinates
            #                 upper_query = DNA(str(altSeq[chrom][Ustart -1 : Ustop - 1 ]))
            #                 lower_alignment = alignment.local_pairwise_align_ssw(lower_query, subject)
            #                 upper_alignment = alignment.local_pairwise_align_ssw(upper_query, subject)
            #                 pdb.set_trace()
            #                 Laln_len = lower_alignment[1]
            #                 Laln_qcoords = lower_alignment[2][0]
            #                 Laln_scoords = lower_alignment[2][1]
            #                 Ualn_len = upper_alignment[1]
            #                 Ualn_qcoords = upper_alignment[2][0]
            #                 Ualn_scoords = upper_alignment[2][1]
            #                 newStart = pstart + (Laln_scoords[0] - Laln_qcoords[0]) + 1
            #                 dist_from_subject_end = len(subject) - Ualn_scoords[1]
            #                 newStop = pstop - dist_from_subject_end + (len(upper_query) - Ualn_qcoords[1])
            #                 assert newStart < newStop
            #                 SP_ChildCoords[gene][key] = [chrom, newStart, newStop]
            #                 SP_ParCoords[par[0]][key] = parInfo
            #             oldGene = gene
            #             Lstart, Ustart = (start, start) # Reinitialize values for new gene
            #             Lstop, Ustop = (stop, stop)
                        # par1:[[par1,geneA,geneB],[par1,geneC,geneD]]
                        # geneC: [[par1,geneC,geneD],[par2,geneC,geneE]]
                        # Alt_Par_Coords[parent][prog]

# def convertSrc_Parents(SP_Dict, alt_annotation, min_exon_size, srcRef, altRef, buff=0):
#     Src_Par_Coords = {}
#     srcSeq = Fasta(srcRef)
#     altSeq = Fasta(altRef)
    # with open(alt_annotation, 'r') as gff:
    #     for i, line in enumerate(gff):
    #         line = line.strip().split("\t")
    #         if line[0][0] == "#": continue
    #         gene = line[8].split("_")[0].split("=")[1]
    #         if gene in SP_Dict and line[2] == "exon": # Loop over exons in annotation
    #             p += 1
    #             chrom = line[0]
    #             start = int(line[3])
    #             stop = int(line[4])
    #             parent = AP_Dict[gene][0][0]
    #             if stop - start >= min_exon_size: 
    #                 n += 1
    #                 if "oldGene" not in locals(): # Check if oldGene has been initialized for this annotation (i.e. is this the first gene)
    #                     oldGene = gene
    #                     Lstart, Ustart = (start, start)
    #                     Lstop, Ustop = (stop, stop)
    #                 elif oldGene == gene: # Still looking at exons from same gene
    #                     if start < Lstart: # is current exon upstream of current most upstream
    #                         Lstart = start
    #                         Lstop = stop
    #                     elif stop > Ustop: #is current exon downstream of current most downstream
    #                         Ustart = start
    #                         Ustop = stop
    #                 else: # We have moved onto a new gene.  Do alignments here
    #                     if n / p < 0.6: print(f"gene {gene} only has {str(n)} out of {str(p)} exons of sufficient size {str(n/p)}") 
    #                     p, n = (0, 0)

    #                     # get seqs from Lstart: Lstop and Ustart:Ustop and align
    #                     lower_query = DNA(srcSeq[chrom][Lstart - 1 : Lstop - 1]) # Pull sequence from reference using current gene coordinates
    #                     upper_query = DNA(srcSeq[chrom][Ustart -1 : Ustop - 1 ])
    #                     lower_alignment = alignment.local_pairwise_align_ssw(lower_query, subject)
    #                     upper_alignment = alignment.local_pairwise_align_ssw(upper_query, subject)
 
    #                     out.write("\t".join([line[0], str(Lstart), str(Lstop), oldGene + ".L", parent, parCoordDict[parent][2], str(parCoordDict[parent][0]), str(parCoordDict[parent][1])]) + "\n")
    #                     out.write("\t".join([line[0], str(Ustart), str(Ustop), oldGene + ".U", parent, parCoordDict[parent][2], str(parCoordDict[parent][0]), str(parCoordDict[parent][1])]) + "\n")
    #                     oldGene = gene
    #                     Lstart, Ustart = (start, start) # Reinitialize values for new gene
    #                     Lstop, Ustop = (stop, stop)
    #                     if split_size != -9: # To split output files or not
    #                         if m % split_size == 0:
    #                             out.close()
    #                             j += 1
    #                             out = open(outs[k] + str(j), 'w')
    # return(Alt_Par_Coords)

# for par in SP_Dict[gene]:
#                     key = ",".join(par)
#                     parInfo = srcCoords[par[0]]
#                     pstart = parInfo[1] - 1 - buff
#                     pstop =  parInfo[2] - 1 + buff
#                     subject = DNA(str(srcSeq[parInfo[0]][pstart : pstop]).replace("N", ""))
#                     # get seqs from Lstart: Lstop and Ustart:Ustop and align
#                     query = DNA(str(altSeq[chrom][start - 1 : stop - 1]).replace("N", "")) # Pull sequence from reference using current gene coordinates
#                     Alignments = lignment.global_pairwise_align_nucleotide(str(query), str(subject), 2, -1, -3, -1)
#                     rc_Alignments = pairwise2.align.globalms(str(query.complement()), str(subject), 2, -1, -3, -1)
#                     sorted_alignments = sorted(Alignments, key=operator.itemgetter(2))
#                     pdb.set_trace()
#                     aln_len = Alignment[1]
#                     rc_aln_len = rc_Alignment[1]
#                     if aln_len / len(query) < threshold and rc_aln_len / len(query) < threshold:
#                         pdb.set_trace()
#                         Alignment = alignment.global_pairwise_align_nucleotide(query, subject)
#                         rc_Alignment = alignment.global_pairwise_align_nucleotide(query.complement(), subject)
#                         pdb.set_trace()
#                     if aln_len < rc_aln_len:
#                         Alignment = rc_Alignment
#                     aln_qcoords = Alignment[2][0]
#                     aln_scoords = Alignment[2][1]
#                     newStart = pstart + (aln_scoords[0] - aln_qcoords[0]) + 1
#                     dist_from_subject_end = len(subject) - aln_scoords[1]
#                     newStop = pstop - dist_from_subject_end + (len(query) - aln_qcoords[1])
                    



            # WILL YOU MISS SOME GENES ENTIRELY 
