'''
    File name: convert_coords_wExons.py
    Author: Patrick Monnahan
    Date created: 09/01/18
    Python Version: 3.6
    Project: Split Genes
    Downstream of: get_boundary_exons.py
    Upstream of: make_psuedo_annotation.py
    Description: For first/last exon pairs in the provided bed file, this script converts coordinates from the original reference to a user-specified reference via automated blasting/parsing. The script get_boundary_exons.py can be used to generate the input
'''

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import argparse
import pdb
import os
import random
import string
import glob


def makeCoordDict(bedfile):
    coord_dict = {} # This dictionary will organize gene/coords within each chromosome
    gene_dict = {} # This dictionary uses genes as keys and simply stores coordinates
    with open(bedfile, 'r') as bed:
        for line in bed:
            line = line.strip("\n").split("\t")
            chrom = line[0].strip("chr0")
            chrom = chrom.strip("chr")
            start = line[1]
            stop = line[2]
            geneID = line[3]
            if geneID.endswith("L"):
                try:
                    gene_dict[geneID[:-2]].append(start)
                except KeyError:
                    gene_dict[geneID[:-2]] = [chrom, start]
            elif geneID.endswith("U"):
                try:
                    gene_dict[geneID[:-2]].append(stop)
                except KeyError:
                    gene_dict[geneID[:-2]] = [chrom, stop]
            try:
                coord_dict[chrom].append([start,stop,geneID])
            except KeyError:
                coord_dict[chrom] = [[start,stop,geneID]]
    # pdb.set_trace()
    return(coord_dict, gene_dict)

def getSequences(coord_dict, source_ref, temp_file):
    with open(temp_file, 'w') as multi_fasta: # We will get sequences corresponding to gene boundaries from the current reference that we wish to blast and store these in a file
        for i, chrom in enumerate(SeqIO.parse(source_ref, "fasta")): # Use BIOpython to parse reference fasta.
            try:
                for coord in coord_dict[chrom.id]: # Loop over all gene coordinates of interest for the current chromosome
                    seq = str(chrom.seq)[int(coord[0]):int(coord[1])] # Pull sequence from reference using current gene coordinates
                    N = str(len(seq) - seq.count("N")) # Count number of N's in the sequence in order to adjust query length that is included in multi_fasta 
                    # pdb.set_trace()
                    if coord[-1].endswith("L"): # L stands for Lower boundary or Lower exon
                        multi_fasta.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2][:-2]}S;{N}\n{seq}\n") # Write tmp gene id and exon sequence
                    else: # Catch genes that end with U for Upper boundary
                        multi_fasta.write(f">{chrom.id}:{coord[0]}-{coord[1]};{coord[2][:-2]}E;{N}\n{seq}\n")
            except KeyError: pass
    multi_fasta.close()
    return()

def Blast(multi_fasta, Toblastdb, cluster=False):
    if cluster: 
        blastn_cline = NcbiblastnCommandline(query=multi_fasta, db=Toblastdb, evalue=0.001,
                                        outfmt=5, max_target_seqs=4, max_hsps_per_subject=4, out=multi_fasta.replace(".fasta",".xml"))
    else:
        blastn_cline = NcbiblastnCommandline(query=multi_fasta, db=Toblastdb, evalue=0.001,
                                        outfmt=5, max_target_seqs=6, max_hsps=6, out=multi_fasta.replace(".fasta",".xml"))
    blastn_cline()
    results = NCBIXML.parse(open(multi_fasta.replace(".fasta",".xml")))
    return(results)

def parseBlastResults(blastResults, newCoords, thresshold, breakFirst=False):
    noHits =[] # Initially stores all results but we will subtract out the found genes to just return those with no hits
    found = []
    for result in blastResults: # each 'result' is for a blasted exon
        noHits.append(result.query)
        gene = result.query.split(";")[1]
        qlen = int(result.query.split(";")[-1]) # exon length - # of N's
        chrom = result.query[0]
        for alignment in result.alignments: # each alignment 
            for hsps in alignment.hsps: # Loop over High Scoring Sequence Pairs with each alignement
                if abs(hsps.align_length - qlen) / qlen < thresshold: # Percent difference of alignment length WRT query must be within thresshold
                    found.append(result.query)
                    if gene[-1] == "S":
                        pos = hsps.sbjct_start - hsps.query_start + 1 # sbjct_start is blast location of ref, but need to modify recorded position based on where first exon ought to start??
                    elif gene[-1] == "E":
                        pos = hsps.sbjct_end + (qlen - hsps.query_end)
                    else: print("Did not find gene bound of query gene") # This should never trip
                    try:
                        newCoords[gene].append(int(pos))
                    except KeyError:
                        newCoords[gene] = [int(chrom), int(pos)]
                    if breakFirst: break # Can be used to just take the first sufficiently matching alignment
            if breakFirst: break

    noHits = list(set(noHits) - set(found)) # Determine which genes were not found by blasting
    return(newCoords, noHits, found)

def removeDups(newCoords, gene_dict):
    """In the case that an exon returns multiple sufficiently good matches (happens frequently), this function will look across all candidate positions for start and end coordinates and select the pair that best matches the expected gene length"""
    dups = []
    for gene in newCoords: # Find duplicates as gene entries with 
        if len(newCoords[gene]) > 2: # items for a gene key include at least the chromosome and one position entry.  Anything greater than two means multiple position entries, hench duplicate
            dups.append(gene)
    # pdb.set_trace()
    done =[] 
    finCoords = {}
    for dup in dups: #Remove duplicate entries based on which entry is closest in expected position
        # pdb.set_trace()
        gene = dup[:-1] # Remove the S or E suffix
        if gene not in done: # Skip duplicate removal if already done for either E or S entry for a particular gene.  
            try:
                starts = newCoords[gene + "S"][1:] # Get list of candidate starts
                ends = newCoords[gene + "E"][1:]
            except KeyError:
                print(f"Did not find new coordinates for {gene}")
                continue # Skip this gene if either E or S are not found
            elen = abs(int(gene_dict[gene][1]) - int(gene_dict[gene][2])) # Expected length of the gene based on old coordinates
            best_match = 999999999 # Initialize for comparison
            # pdb.set_trace()
            for i, start in enumerate(starts):
                olens = [abs(start - end) for end in ends] # Calculate the observed length for every candidate end with current start
                len_diffs = [abs(olen - elen) for olen in olens]
                new_best = min(len_diffs) # Find end position that minimizes difference in observed and expected length for current start
                if new_best < best_match: # See if current start has a candidate end that provides a better match to expected length
                    best_match = new_best
                    idx = len_diffs.index(new_best)
                    winners = [i + 1, idx + 1] # Record index positions within start and end lists for the selected best match coordinates
            # pdb.set_trace()
            finCoords[gene] = [newCoords[gene + "S"][0], newCoords[gene + "S"][winners[0]], newCoords[gene + "E"][winners[1]]] # Report the identified winners from the S and E lists after going through all possible start and end pairs    
            done.append(gene)
    for gene in newCoords: # Complete final coordinate list by adding all non-duplicated genes
        gene = gene[:-1]
        if gene not in done:
            try:
                finCoords[gene] = [newCoords[gene + "S"][0], newCoords[gene + "S"][1], newCoords[gene + "E"][1]]
            except KeyError:
                print(f"Did not find new coordinates for {gene}")
    return(finCoords)

def getNewCoords(multi_fasta, Toblastdb, thresshold, gene_dict, cluster=False):
    # pdb.set_trace()
    results = Blast(multi_fasta, Toblastdb, cluster)
    newCoords ={}
    newCoords, noHits, found = parseBlastResults(results, newCoords, thresshold)
    # pdb.set_trace()
    newCoords = removeDups(newCoords, gene_dict)
    # pdb.set_trace()
    return(newCoords)

def verifyNewCoords(newCoords, gene_dict, threshold):
    """This function checks whether the new gene lengths determined via blasting are within a given threshold of the old distance determined from the annotation"""
    verCoords = {}
    for gene, coords in gene_dict.items():
        # pdb.set_trace()
        start = int(coords[1])
        stop = int(coords[2])
        length = abs(stop - start)  ##THIS HAS TO ABS IN ORDER TO ACCOUNT FOR ALIGNMENT TO OPPOSITE STRAND
        try:
            nchrom, nstart, nstop = newCoords[gene] # info from new coordinates
            nlen = abs(nstop - nstart)
            if abs(nlen - length) / length > threshold: # Difference in new and old length expressed as a percent of the old length
                print(f"bad match for {gene}: old {start}-{stop} {int(stop)-int(start)} new {nstart}-{nstop} {nlen}")
            else:
                if stop - start < 0: # Write entry such that start precedes stop
                    verCoords[gene] = [str(newCoords[gene][0]), str(newCoords[gene][2]), str(newCoords[gene][1])]
                else:
                    verCoords[gene] = [str(x) for x in newCoords[gene]]
                print(f"good match for {gene}: old {start}-{stop} {int(stop)-int(start)} new {nstart}-{nstop} {int(nstop) - int(nstart)}")
        except KeyError: pass
    return(verCoords)

def writeOutBed(outFile, newCoords, oldCoords):
    """newCoords should be output of verifyNewCoords and oldCoords should be gene_dict from makeCoordDict"""
    with open(outFile, 'w') as out: # Write output to file
        for gene in verCoords:
            new_coord = "\t".join(newCoords[gene])
            old_coord = "\t".join(oldCoords[gene])
            new_len = int(verCoords[gene][2]) - int(verCoords[gene][1])
            old_len = int(gene_dict[gene][2]) - int(gene_dict[gene][1])
            out.write(f"{new_coord}\t{gene}\t{old_coord}\t{new_len}\t{old_len}\n")
    return()

if __name__ == "__main__":
    import time
    start_time = time.time()
    parser = argparse.ArgumentParser(description = "For first/last exon pairs in the provided bed file, this script converts coordinates from the original reference to a user-specified reference via automated blasting/parsing. The script get_boundary_exons.py can be used to generate the input")
    parser.add_argument('-b', type=str, metavar='bed_file', required=True, help='bed file containing coordinates and geneIDs for genes you wish to convert')
    parser.add_argument('-d', type=str, metavar='blastdb', required=True, help='Full path to blast database.  This corresponds to the reference whose coordinates you are trying to convert to')
    parser.add_argument('-r', type=str, metavar='source_ref', required=True, help='Reference fasta corresponding to reference converting from as in bed file')
    parser.add_argument('-t', type=str, metavar="tmp_dir", default=os.getcwd(), help="temporary directory that will store multi fasta for blasting along with blast results")
    parser.add_argument('-T1', type=int, metavar="size_similarity_threshold", default=0.5, help="percent difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the query length")
    parser.add_argument('-T2', type=int, metavar="size_similarity_threshold", default=3.0, help="proportion difference between gene size as function of old and new coordinates")
    parser.add_argument('-O', type=str, metavar='output_bed_name', required=True, help="bed file to write output to.  Fields are new chrom, start and stop position followed by geneID, the old coordinates, and the new and old gene lengths")
    parser.add_argument('-v', action="store_true")
    parser.add_argument('--cluster', action="store_true", help="set this flag if running on cluster.  Necessary for using MSI blast installations")
    parser.add_argument('-k', action="store_true", help="keep temporary files: fastas and blast results")
    args = parser.parse_args()


    coord_dict, gene_dict = makeCoordDict(args.s) # Parse info in input bed file
    print("Made Coordinate Dictionary")
    tmp_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)) # Make random string for tmp fasta and blast results
    tmp_file = f"{args.t}/tmp{tmp_str}.fasta" 
    getSequences(coord_dict, args.r, tmp_file) # Retrieve original exon sequences
    print("Retrieved sequences of coordinates")
    newCoords = getNewCoords(tmp_file, args.b, args.T1, gene_dict, args.cluster) # get new coordinates via blasting exon sequence
    # pdb.set_trace()
    verCoords = verifyNewCoords(newCoords, gene_dict, args.T2) # Verify that new gene lengths are within given distance from old gene lengths
    writeOutBed(args.O, verCoords, gene_dict) # Write output to bed file

    if not args.k:
        for j in glob.glob(f'{args.t}/tmp{tmp_str}*'):
            os.remove(j)

