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

from Bio import SeqIO
import argparse
import pdb
from skbio import DNA, alignment

def makeCoordDict(bedfile):
    lower_exon_dict = {str(x): {} for x in range(1,11)} # This dictionary uses chromosomes as keys and stores each exon coordinates along with parent gene + coords
    upper_exon_dict = {str(x): {} for x in range(1,11)}
    with open(bedfile, 'r') as bed:
        for line in bed:
            chrom, start, stop, geneID, parID, pchrom, pstart, pstop = line.strip("\n").split("\t")
            Chrom = chrom.lstrip("chr0").strip("chr")
            Pchrom = pchrom.lstrip("chr0").strip("chr")
            assert Chrom == Pchrom
            if geneID.endswith("L"): 
                lower_exon_dict[Chrom][geneID[:-2]] = [chrom, start, stop, geneID, pchrom, pstart, pstop, parID]
            else:
                upper_exon_dict[Chrom][geneID[:-2]] = [chrom, start, stop, geneID, pchrom, pstart, pstop, parID]
    # pdb.set_trace()
    return(lower_exon_dict, upper_exon_dict)

def convertCoords(lower_exon_dict, upper_exon_dict, srcRef, altRef, buff=0):
    newCoords = {}
    for i, chrom in enumerate(SeqIO.parse(srcRef, "fasta")): # Use BIOpython to parse reference fasta.
        for j, pchrom in enumerate(SeqIO.parse(altRef, "fasta")):
            Chrom = chrom.id.strip("chr0").strip("chr")
            Pchrom = pchrom.id.strip("chr0").strip("chr")
            if Chrom == Pchrom:
                Alt = pchrom
        for gene, Linfo in lower_exon_dict[chrom.id].items(): # Loop over all gene coordinates of interest for the current chromosome
            Uinfo = upper_exon_dict[chrom.id][gene]
            pdb.set_trace()
            # c, lstart, lstop, geneID, p, pstart, pstop, laltID = lower_exon_dict[chrom.id][gene]
            # c, ustart, ustop, geneID, p, pstart, pstop, ualtID = 
            lower_query = DNA(str(chrom.seq)[ int(Linfo[1]) : int(Linfo[2]) ]) # Pull sequence from reference using current gene coordinates
            upper_query = DNA(str(chrom.seq)[ int(Uinfo[1]) : int(Uinfo[2]) ]) # Pull sequence from reference using current gene coordinates
            assert Linfo[-1] == Uinfo[-1]
            subject = str(Alt.seq)[ int(Linfo[5]) - buff : int(Linfo[6]) + buff ] # Pull sequence from reference using current gene coordinates
            sN = str(len(subject) - subject.count("N")) # Count number of N's in the sequence in order to adjust query length that is included in multi_fasta 
            lower_alignment = alignment.local_pairwise_align_ssw(lower_query, subject)
            upper_alignment = alignment.local_pairwise_align_ssw(upper_query, subject)
            Laln_len = lower_alignment[1]
            Laln_qcoords = lower_alignment[2][0]
            Laln_scoords = lower_alignment[2][1]
            Ualn_len = upper_alignment[1]
            Ualn_qcoords = upper_alignment[2][0]
            Ualn_scoords = upper_alignment[2][1]

            propL_aln = abs(Lstop - Lstart) / len(lower_query)
            propU_aln = abs(Ustop - Ustart) / len(upper_query)

            ## STOPPED HERE...DEAL with parsing alignment results
            ## HOW TO CONVERT ALIGNMENT RESULTS TO NEW COORDINATES?? 
            ## DRAW OUT ALL POSSIBLE WAYS IN WHICH THEY COULD ALIGN!!!
            pdb.set_trace()
            if Ustart > Lstart: # Same orientation between src and alt
                assert Lstart < Lstop and Ustart < Ustop
                newStart = int(pstart) - buff + Lstart
                newStop = int(pstop) + buff - Ustop
            else:
                assert Lstart > Lstop and Ustart > Ustop
                newStop = int(pstart) - buff + Ustart
                newStart = int(pstop) + buff - Lstart
            oldStart = min(lstart, ustart)
            oldStop = max(lstop, ustop)
            newCoords[gene] = [p, newStart, newStop, geneID, (newStop - newStart), c, oldStart, oldStop, oldStop - oldStart]
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

def writeOutBed(outFile, newCoords):
    """newCoords should be output of verifyNewCoords and oldCoords should be gene_dict from makeCoordDict"""
    with open(outFile, 'w') as out: # Write output to file
        for gene in newCoords:
            new_coord = "\t".join(newCoords[gene]) + "\n"
            out.write(new_coord)
    return()

if __name__ == "__main__":
    import time
    start_time = time.time()
    parser = argparse.ArgumentParser(description = "For first/last exon pairs in the provided bed file, this script converts coordinates from the original reference to a user-specified reference via automated blasting/parsing. The script get_boundary_exons.py can be used to generate the input")
    parser.add_argument('-b', type=str, metavar='bed_file', required=True, help='bed file containing coordinates and geneIDs for genes you wish to convert')
    parser.add_argument('-d', type=str, metavar='blastdb', required=True, help='Full path to alternative reference whose coordinates you are trying to convert to.  If using --global, then you must run makeblastdb on this fasta prior to running')
    parser.add_argument('-r', type=str, metavar='source_ref', required=True, help='Reference fasta corresponding to reference converting FROM as in bed file')
    parser.add_argument('-T1', type=int, metavar="size_similarity_threshold", default=0.5, help="percent difference between aligned length and query length must be within this thresshold, which is specified as a proportion of the query length")
    parser.add_argument('-T2', type=int, metavar="size_similarity_threshold", default=3.0, help="proportion difference between gene size as function of old and new coordinates")
    parser.add_argument('-B', type=int, metavar="buffer", default=2000, help="For local blast to expected gene, how much buffer to add to both sides of the sbjct sequence")
    parser.add_argument('-O', type=str, metavar='output_bed_name', required=True, help="bed file to write output to.  Fields are new chrom, start and stop position followed by geneID, the old coordinates, and the new and old gene lengths")
    parser.add_argument('-v', action="store_true")
    args = parser.parse_args()

    LED, UED = makeCoordDict(args.b) # Parse info in input bed file for parent gene name in 5th column
    print("Made Coordinate Dictionary for alternative reference")
    new = convertCoords(LED, UED, args.r, args.d, args.B)
    print("Retrieved sequences of coordinates from alternative reference")

    # pdb.set_trace()
    # verCoords = verifyNewCoords(newCoords, gene_dict, args.T2) # Verify that new gene lengths are within given distance from old gene lengths
    writeOutBed(args.O, new, gene_dict) # Write output to bed file

