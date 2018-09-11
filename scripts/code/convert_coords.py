from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import argparse
import pdb
import os
import random
import string
import collections
from math import ceil
import bx.align.lav
import glob

# a = bx.align.lav.LavAsPiecesReader(open("/Users/pmonnahan/Documents/Research/MaizeSV/junk/testlast.out.lav",'rb'))

# def Lastz(geneID, multi_fasta, refpath, lastz_path):
#     cmd = f"grep -A1 {geneID} {multi_fasta} | {lastz_path} {refpath}"


def makeCoordDict(bedfile):
    coord_dict = {}
    gene_dict = {}
    with open(bedfile, 'r') as bed:
        for line in bed:
            line = line.strip("\n").split("\t")
            chrom = line[0]
            start = line[1]
            stop = line[2]
            geneID = line[3]
            gene_dict[geneID] = [chrom,start,stop]
            try:
                coord_dict[chrom].append([start,stop,geneID])
            except KeyError:
                coord_dict[chrom] = [[start,stop,geneID]]
    # pdb.set_trace()
    return(coord_dict, gene_dict)


def getSequences(coord_dict, source_ref, temp_file, Offset):
    fullSequences = open(temp_file.replace(".fasta",".full.fasta"), 'w')
    with open(temp_file, 'w') as multi_fasta:
        for i, chrom in enumerate(SeqIO.parse(source_ref, "fasta")):
            try:
                for coord in coord_dict[chrom.id]:
                    if int(coord[1]) - int(coord[0]) < Offset:
                        offset = ceil(abs(0.2 * (int(coord[0]) - int(coord[1]))))
                    else: offset = Offset
                    seq1 = str(chrom.seq)[int(coord[0]):int(coord[0]) + offset]
                    seq2 = str(chrom.seq)[int(coord[1]) - offset:int(coord[1])]
                    Nstart = str(len(seq1) - seq1.count("N"))
                    Nstop = str(len(seq2) - seq2.count("N"))
                    multi_fasta.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2]}S;{Nstart}\n{seq1}\n")
                    multi_fasta.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2]}E;{Nstop}\n{seq2}\n")
                    full = str(chrom.seq)[int(coord[0]):int(coord[1])]
                    fullSequences.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2]}F;{len(full)}\n{full}\n")
            except KeyError: pass
    multi_fasta.close()
    return()

def getBackupSeqs(fullSeqFasta, queryIDlist, numbkps = 5):
    gene_dict = {}
    for query in queryIDlist:
        gene = query.split(";")[1]
        chrom = query[0]
        loc = gene[-1]
        offset = int(query.split(";")[-1])
        gene_dict[gene[:-1]] = [chrom, loc, offset]
        
    bkp_fasta = open(fullSeqFasta.replace(".full.fasta", ".bkp.fasta"), 'w')
    for i, entry in enumerate(SeqIO.parse(fullSeqFasta, "fasta")):
        # pdb.set_trace()
        gene = entry.id.split(";")[1][:-1]
        if gene in gene_dict or gene in gene_dict:
            chrom, loc, Offset = gene_dict[gene]
            if loc == 'E':
                pos2 = len(entry.seq)
                pos1 = pos2 - Offset
                for i in range(0,numbkps):
                    pos2 -= Offset
                    pos1 -= Offset
                    Seq = str(entry.seq)[pos1:pos2]
                    if pos2 > 0 and pos1 > 0:
                        qlen = str(len(Seq) - Seq.count("N"))
                        bkp_fasta.write(f">{chrom}:{pos1}-{pos2};{gene}E;{qlen}\n{Seq}\n")
            elif loc == 'S':
                pos2 = Offset
                pos1 = 0
                for i in range(0,numbkps):
                    pos2 += Offset
                    pos1 += Offset
                    Seq = str(entry.seq)[pos1:pos2]
                    if Seq:
                        qlen = str(len(Seq) - Seq.count("N"))
                        bkp_fasta.write(f">{chrom}:{pos1}-{pos2};{gene}S;{qlen}\n{Seq}\n")
    bkp_fasta.close()
    return(fullSeqFasta.replace(".full.fasta", ".bkp.fasta"))

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

def Blast(multi_fasta, Toblastdb, cluster=False):
    if cluster:
        blastn_cline = NcbiblastnCommandline(query=multi_fasta, db=Toblastdb, evalue=0.001,
                                        outfmt=5, max_target_seqs=4, max_hsps_per_subject=4, out=multi_fasta.replace(".fasta",".xml"))
    else:
        blastn_cline = NcbiblastnCommandline(query=multi_fasta, db=Toblastdb, evalue=0.001,
                                        outfmt=5, max_target_seqs=4, max_hsps=4, out=multi_fasta.replace(".fasta",".xml"))
    blastn_cline()
    results = NCBIXML.parse(open(multi_fasta.replace(".fasta",".xml")))
    return(results)

def removeDups(newCoords, gene_dict):
    dups = []
    for gene in newCoords:
        if len(newCoords[gene]) > 2:
            dups.append(gene)
    # pdb.set_trace()
    for dup in dups: #Remove duplicate entries based on which entry is closest in expected position
        # pdb.set_trace()
        if dup[-1] == "S":
            comp = gene_dict[dup[:-1]][1]
        else:
            comp = gene_dict[dup[:-1]][2]
        values = [abs(d - int(comp)) for d in newCoords[dup][1:]]
        idx = values.index(min(values)) + 1
        newCoords[dup] = [newCoords[dup][0], newCoords[dup][idx]]
    return(newCoords)


def getNewCoords(multi_fasta, Toblastdb, thresshold, gene_dict, thresh2, thresh3, cluster=False):
    
    # pdb.set_trace()
    results = Blast(multi_fasta, Toblastdb, cluster)
    newCoords ={}
    newCoords, noHits, found = parseBlastResults(results, newCoords, thresshold)
    # pdb.set_trace()
    if noHits: #Is this right way to deal with no hits?
        bkp_fasta = getBackupSeqs(multi_fasta.replace(".fasta",".full.fasta"), noHits)
        newResults = Blast(bkp_fasta, Toblastdb, cluster)
        newCoords, noHits2, found = parseBlastResults(newResults, newCoords, thresh2, breakFirst=True)
    # pdb.set_trace()

    fo = [query.split(";")[1] for query in found]
    sec = [query.split(";")[1] for query in noHits2]
    if list(set(sec) - set(fo)):
        print(f"Did not find hits for {sec}")
    # pdb.set_trace()

    newCoords = removeDups(newCoords, gene_dict)
    finCoords = {}
    for ID in newCoords:
        # pdb.set_trace()
        gene = ID[:-1]
        try:
            finCoords[gene] = [newCoords[ID][0], newCoords[gene + "S"][1], newCoords[gene + "E"][1]]
        except KeyError:
            print(f"Did not find new coordinates for {gene}")
    return(finCoords)

def verifyNewCoords(newCoords, bed_file, threshold):
    verCoords = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip("\n").split("\t")
            # pdb.set_trace()
            start = int(line[1])
            stop = int(line[2])
            length = abs(stop - start)  ##THIS HAS TO ABS IN ORDER TO ACCOUNT FOR ALIGNMENT TO OPPOSITE STRAND
            geneID = line[3] 
            try:
                nchrom, nstart, nstop = newCoords[geneID]
                nlen = abs(nstop - nstart)
                if abs(nlen - length) / length > threshold:
                    print(f"bad match for {geneID}: old {start}-{stop} {int(stop)-int(start)} new {nstart}-{nstop} {nlen}")
                else:
                    if stop - start < 0:
                        verCoords[geneID] = [str(newCoords[geneID][0]), str(newCoords[geneID][2]), str(newCoords[geneID][1])]
                    else:
                        verCoords[geneID] = [str(x) for x in newCoords[geneID]]
                    print(f"good match for {geneID}: old {start}-{stop} {int(stop)-int(start)} new {nstart}-{nstop} {int(nstop) - int(nstart)}")
            except KeyError: pass
    return(verCoords)


if __name__ == "__main__":
    import time
    start_time = time.time()
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument('-s', type=str, metavar='bed_file_of_genes_to_convert', required=True, help='')
    parser.add_argument('-b', type=str, metavar='blastdb', required=True, help='Full path to blast databases')
    parser.add_argument('-r', type=str, metavar='source_ref', required=False, default = "/Users/pmonnahan/Documents/Research/MaizeSV/references/B73_genesOnly.bed", help='Full path to GSTRiP_RedundancyAnnotator.sh')
    parser.add_argument('-t', type=str, metavar="tmp_dir", default=os.getcwd())
    parser.add_argument('-o', type=int, metavar="offset", default=1000, help="This determines the length of the sequence from the start and end positions to be used for blast to get alt ref coordinates")
    parser.add_argument('-T1', type=int, metavar="size_similarity_threshold", default=0.5, help="difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the offset value")
    parser.add_argument('-T2', type=int, metavar="size_similarity_threshold", default=0.8, help="difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the offset value")
    parser.add_argument('-T3', type=int, metavar="size_similarity_threshold", default=0.7, help="difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the offset value")
    parser.add_argument('-T4', type=int, metavar="max_distance_to_alt_pos", default=2000000, help="difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the offset value")
    parser.add_argument('-O', type=str, metavar='output_bed_name', required=True)
    parser.add_argument('-v', action="store_true")
    parser.add_argument('--cluster', action="store_true", help="set this flag if running on cluster.  Necessary for using MSI blast installations")
    parser.add_argument('-k', action="store_true", help="keep temporary files: fastas and blast results")
    args = parser.parse_args()


    coord_dict, gene_dict = makeCoordDict(args.s)
    print("Made Coordinate Dictionary")
    tmp_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    tmp_file = f"{args.t}/tmp{tmp_str}.fasta" 
    getSequences(coord_dict, args.r, tmp_file, args.o)
    print("Retrieved sequences of coordinates")
    newCoords = getNewCoords(tmp_file, args.b, args.T1, gene_dict, args.T2, args.T4, args.cluster)
    verCoords = verifyNewCoords(newCoords, args.s, args.T3)
    # pdb.set_trace()
    with open(args.O, 'w') as out:
        for gene in verCoords:
            new_coord = "\t".join(verCoords[gene])
            old_coord = "\t".join(gene_dict[gene])
            new_len = int(verCoords[gene][2]) - int(verCoords[gene][1])
            old_len = int(gene_dict[gene][2]) - int(gene_dict[gene][1])
            out.write(f"{new_coord}\t{gene}\t{old_coord}\t{new_len}\t{old_len}\n")
    if not args.k:
        for j in glob.glob(f'{args.t}/tmp{tmp_str}*'):
            os.remove(j)


    # print("--- %s seconds ---" % (time.time() - start_time))



        # unfound = set(genes) - set(found)
    # # pdb.set_trace()
    # attempts = 0
    # while unfound:
    #     attempts += 1
    #     if attempts > 10:
    #         print(f"Unable to find matches for {unfound}")
    #         pdb.set_trace()
    #         break
    #     results = NCBIXML.parse(open(multi_fasta.replace(".fasta",".xml")))
    #     for result in results:
    #         gene = result.query.split(";")[1]
    #         if gene in unfound:
    #             found = False
    #             if gene[:-1] + "E" in newCoords or gene[:-1] + "S" in newCoords: #was partner found that we can use as an anchor
    #                 try:
    #                     anchor = newCoords[gene[:-1] + "E"][1]
    #                 except KeyError:
    #                     anchor = newCoords[gene[:-1] + "S"][1]
    #                 # pdb.set_trace()
    #                 eCoords = gene_dict[gene[:-1]]
    #                 eDist = abs(int(eCoords[1]) - int(eCoords[2]))
    #                 qlen = int(result.query.split(";")[-1])
    #                 chrom = result.query[0]
    #                 for alignment in result.alignments: 
    #                     for hsps in alignment.hsps:
    #                         # pdb.set_trace()
    #                           # Stopped here...how to find right value in newCoords to calculate distance relative to???
    #                         dist = max([abs(hsps.sbjct_start - anchor), abs(hsps.sbjct_end - anchor)])
    #                         if (abs((dist - eDist)) / eDist) < thresh2:
    #                             if gene[-1] == "S":
    #                                 pos = hsps.sbjct_start - hsps.query_start + 1 #plus one because query_start will start at 1
    #                             elif gene[-1] == "E":
    #                                 pos = hsps.sbjct_end + (qlen - hsps.query_end)
    #                             newCoords[gene] = [chrom, pos]
    #                             unfound.remove(gene)
    #                             found = True
    #                             break
    #                     if found: break
    #             else: print(f"Did not find anchoring partner for {gene}")                
    #     print(f"did not find good matches for: {','.join(unfound)}")
    # Sort the entries in the dictionary prior to returning, such that the order will be chr, start, stop
        # if dup[-1] == "S":
        #     anchor = np.mean(newCoords[dup[:-1] + "E"][1:])

        # else:
        #     anchor = np.mean(newCoords[dup[:-1] + "S"][1:])
        # dist = [abs(anchor - d) for d in newCoords[dup][1:]]
        # dist2 = [abs(d - eDist) for d in dist]
        # idx = dist2.index(min(dist2)) + 1
        # newCoords[dup] = [newCoords[dup][0], newCoords[dup][idx]]

    #                 for alignment in result.alignments: 
    #                     for hsps in alignment.hsps:
    #                         # pdb.set_trace()
    #                           # Stopped here...how to find right value in newCoords to calculate distance relative to???
    #                         dist = max([abs(hsps.sbjct_start - anchor), abs(hsps.sbjct_end - anchor)])
    #                         if (abs((dist - eDist)) / eDist) < thresh2: 
        #     comp = gene_dict[dup[:-1]][1]
        # else:
        #     comp = gene_dict[dup[:-1]][2]
        # values = [abs(d - int(comp)) for d in newCoords[dup][1:]]
        # idx = values.index(min(values)) + 1
        # newCoords[dup] = [newCoords[dup][0], newCoords[dup][idx]]
                # eCoords = gene_dict[dup[:-1]]
        # eDist = abs(int(eCoords[1]) - int(eCoords[2]))
