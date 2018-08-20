import pickle
import numpy as np
import allel
import argparse
import pdb
import line_profiler
import subprocess
import csv


# @profile
def makeBedDict(file_list, verbose=False):
    coord_dict = {}
    for bed_file in file_list:
        with open(bed_file, 'r') as bed:
            for i, line in enumerate(bed):
                line = line.strip("\n").split("\t")
                # pdb.set_trace()
                if i == 0:
                    coord_dict[line[3][6]] = {}
                # Chromosome names must only be numbers and not include 'chr'
                try:
                    coord_dict[line[3][6]][line[3]] = line[0:3] # line[3][6] finds the reference identifier within the gene name
                except KeyError:
                    if verbose: print(f"bad keys for coord_dict with either {line[3][6]} or {line[3]}")
    # pdb.set_trace()
    return coord_dict

def makeSplitGeneKeyDict(gene_key_file):
    with open(gene_key_file, 'r') as split_gene_key:
        parent_list = []
        SGK = {}
        for i, line1 in enumerate(csv.reader(split_gene_key, delimiter = ",")):
            parent = line1[0]
            children = sorted(list(filter(None, line1[1:])))
            cref = children[0][6]
            if i == 0:
                SGK[parent] = [children]
            elif parent in parent_list:
                if SGK[parent] == children:
                    if args.v: print(f"Duplicate entry for {SGK[parent]} and {children} parent: {parent}")
                elif SGK[parent][0][6] == cref:
                    if args.v: print(f"Multiple sets of children from same reference for {SGK[parent]} and {children} parent: {parent}")
                else:
                    SGK[parent].append(children)
            else:
                SGK[parent] = children
            parent_list.append(parent)
    return(SGK)

def newGenes(ref1name, prog_set, coord_dict):
    MetaKey = {"1": "b", "8": "p", "4": "w"}
    pref = parent[6]
    pid = parent[7:]
    childDict = {}
    metaInfo = [MetaKey[pref].upper() + "1"]
    for ref2name_list in prog_set:
        cref = ref2name_list[0][6]
        cnum = len(ref2name_list)
        metaInfo.append(MetaKey[cref] + str(cnum)) 
    if len(metaInfo) != 3:
        oref = [k for k in ["b","p","w"] if k not in [MetaKey[pref], MetaKey[cref]]]
        metaInfo.append(oref[0] + "0")
    metaInfo = "".join(sorted(metaInfo, key=str.lower))
    for ref2name_list in prog_set:
        for child in ref2name_list:
            newID = metaInfo + pid + child[7:]
            childDict[newID] = coord_dict[cref][child]
    pID = metaInfo + pid + pid.upper()   
    parDict = {pID: coord_dict[pref][parent]}
    return(parDict, childDict)

def updateGeneID(ref1name, prog_set, coord_dict, oldname):
    MetaKey = {"1": "b", "8": "p", "4": "w"}
    pref = parent[6]
    pid = parent[7:]
    childDict = {}
    metaInfo = [MetaKey[pref].upper() + "1"]
    for ref2name_list in prog_set:
        cref = ref2name_list[0][6]
        cnum = len(ref2name_list)
        metaInfo = metaInfo.append(MetaKey[cref] + str(cnum)) 
    if len(metaInfo) != 3:
        oref = [k for k in ["b","p","w"] if k not in [MetaKey[pref], MetaKey[cref]]]
        metaInfo.append(oref[0] + "0")
    metaInfo = sorted(metaInfo, key=str.lower)
    for ref2name_list in prog_set:
        for child in ref2name_list:
            newID = metaInfo + pid + child[7:]
            childDict[newID] = coord_dict[cref][child]
    pID = metaInfo + pid + pid.upper()   
    parDict = {pID: coord_dict[pref][parent]}
    return(parDict, childDict)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument('-s', type=str, metavar='split-gene-key', required=True, help='comma separated file where first entry is the parent/merged gene and all subsequent entries are corresponding genes from split model in alternative reference')
    parser.add_argument('-b', type=str, metavar='b73_bed', required=False, default = "/Users/pmonnahan/Documents/Research/MaizeSV/references/B73_genesOnly.bed", help='Full path to GSTRiP_RedundancyAnnotator.sh')
    parser.add_argument('-w', type=str, metavar='w22_bed', required=False, default = "/Users/pmonnahan/Documents/Research/MaizeSV/references/W22_genesOnly.bed")
    parser.add_argument('-p', type=str, metavar='ph207_bed', required=True, help = "/Users/pmonnahan/Documents/Research/MaizeSV/references/PH207_genesOnly.bed")
    parser.add_argument('-v', action="store_true")
    # parser.add_argument('-o', type=str, metavar='output_VCF_name', required=True, help = "Full path to final merged vcf")
    # parser.add_argument('-k', action='store_true', help="Keep temporary files")
    args = parser.parse_args()

    coord_dict = makeBedDict([args.b, args.w, args.p], args.v)
    # print(coord_dict)
    
    SGK = makeSplitGeneKeyDict(args.s)

    mergeds = {}
    splits = {}
    parent_list = []
    child_list = []

    for parent in SGK:
        Nested = False
        for prog_set in SGK[parent]:
            if any(k in parent_list for k in prog_set):
                # Rename grandparent
                print(f"Nested event found for child/parent {[k for k in prog_set if k in parent_list]} and (grand)parent {parent}")
                Nested = True
                # What about imperfect nesting events where grandchildren do not all belong to grandparent??
            if any(k in child_list for k in prog_set):
                # Rename grandparent
                print(f"Child gene(s) {[k for k in prog_set if k in child_list]} has multiple parents in addition to {parent}")
        if not Nested:
            print(SGK[parent])
            parGene, childGenes = newGenes(parent, SGK[parent], coord_dict)
            splits = dict(splits, **childGenes)
            mergeds = dict(mergeds, **parGene)
            print(splits)
            for prog_set in SGK[parent]:
                child_list = child_list + prog_set
        parent_list.append(parent)
        

    # MUST ALSO DEAL WITH DUPLICATES WHOSE CHILDREN MAY BE SORTED DIFFERENTLY ACROSS DUPLICATES
