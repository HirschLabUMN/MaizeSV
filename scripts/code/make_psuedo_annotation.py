'''
    File name: make_psuedo_annotation.py
    Author: Patrick Monnahan
    Date created: 09/01/18
    Python Version: 3.6
    Project: Split Genes
    Downstream of: convert_coords_wExons.py
    Upstream of: ??
    Description: This script takes a series of split/merged candidates and produces a single pseudo-annotation file containing coordinates of each gene along with their corresponding split/merged genes.  The coordinates within the bed files (input in -b) should all be in terms of a common reference as accomplished with convert_coords_wExons.py
'''

import argparse
import pdb
import csv

MetaKey = {"1": "b", "8": "p", "4": "w"} # Used to link gene names (e.g. Zm00001d00174) to reference letter: 1 = b73, 8 = ph207, and 4 = w22
MetaKey2 = {"d": 0, "a": 2, "b": 4} # Used to link the characters letters in gene ID suffixes to their index position in the metaInfo string for new gene naming
# @profile

def capID(name,x): return (name[0:x] + name[x].upper() + name[x+1:]) # Capitalize xth value in string 'name'
def lowID(name,x): return (name[0:x] + name[x].lower() + name[x+1:]) # xth value in string 'name' will be made lowercase

def find_path(graph, start, end, path=[]):
    # Not currently being used, but is a way to find members of a directed acyclic graph
    path = path + [start]
    if start == end:
        return path
    if not graph.has_key(start):
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None

def makeBedDict(file_list, verbose=False):
    coord_dict = {} # Dictionary of dictionaries with reference ID (e.g 1,4, or 8) as first key and geneID as second level of keys
    for bed_file in file_list:
        with open(bed_file, 'r') as bed:
            for i, line in enumerate(bed):
                line = line.strip("\n").split("\t")
                if i == 0: # New reference bed file.  Must initialize upper level of dictionary
                    coord_dict[line[3][6]] = {}
                # Chromosome names must only be numbers and not include 'chr'
                try:
                    coord_dict[line[3][6]][line[3]] = line[0:3] # line[3][6] finds the reference identifier within the gene name. line[0:3] is the chrom, start, stop and geneID, respectively
                except KeyError: # Did not find the gene ID within the coordinate dictionary.  Means that prior step (i.e convert_coords_wExons.py) did not work for this gene, likely because blast did not find a good match for first or last exon or both.
                    if verbose: print(f"bad keys for coord_dict with either {line[3][6]} or {line[3]}")
    return(coord_dict)

def makeSplitGeneKeyDict(gene_key_file):
    with open(gene_key_file, 'r') as split_gene_key:
        parent_list = []
        SGK = {} # This will be used if creating new gene names via --newNames
        SGK2 = {} # Temporary dictionary for creating bedfiles with traditional names
        for i, line1 in enumerate(csv.reader(split_gene_key, delimiter = ",")):
            parent = line1[0]
            children = sorted(list(filter(None, line1[1:]))) # Remove empty entries in the child entries and sort entries (helps for identifying duplicate entries)
            cref = children[0][6] # Reference ID of child genes.
            try: # Only works if parent gene is already present in dictionary
                SGK2[parent].append([parent] + children)
            except KeyError:
                SGK2[parent] = [[parent] + children]
            for child in children:
                try:
                    SGK2[child].append([parent] + children)
                except KeyError:
                    SGK2[child] = [[parent] + children]
            if i == 0:
                SGK[parent] = [children]
            elif parent in parent_list:
                if any(k == children for k in SGK[parent]):
                    if args.v: print(f"Duplicate entry for {SGK[parent]} and {children} parent: {parent}")
                elif [prog for prog in SGK[parent] if prog[0][6] == cref]:
                    if args.v: print(f"Multiple sets of children from same reference for {SGK[parent]} and {children} parent: {parent}")
                else:
                    SGK[parent].append(children)
            else:
                SGK[parent] = [children]
            parent_list.append(parent)
    SGK3 = {}
    for key in SGK2:
        for i, val in enumerate(SGK2[key]):
            if i == 0:
                SGK3[key] = [val]
            elif val not in SGK3[key]:
                SGK3[key].append(val)
    return(SGK, SGK3)

def newGenes(parent, prog_set, coord_dict, to_coord):
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
    forParCoords = [] 
    for ref2name_list in prog_set: # Get progeny names that match to_coords to pass as key to coord_dict to get parCoords
        for child in ref2name_list:

            cref = ref2name_list[0][6]
            newID = metaInfo + pid + child[7:] + 'c' + to_coord
            try:
                childDict[newID] = coord_dict[cref][lowID(child,7)]
                if cref == to_coord:
                    forParCoords.append(child)
            except KeyError: 
                childDict[newID] = [cref, -9 , -9]
                print(f"Did not find converted coordinates for child {child}")
    pID = metaInfo + pid + pid + 'c' + to_coord
    coord_list = []
    for gene in forParCoords:
        coord_list += [int(coord_dict[to_coord][lowID(gene,7)][1]), int(coord_dict[to_coord][lowID(gene,7)][2])]
    if coord_list:
        parDict = {pID: [coord_dict[to_coord][lowID(forParCoords[0],7)][0], min(coord_list), max(coord_list)]}
    else:
        try:
            parDict = {pID: coord_dict[pref][parent]}
        except KeyError:
            print(f"Did not find converted coordinates for parent {parent}")
            parDict = {pID: [-9,-9,-9]}
    return(parDict, childDict)

def modDict(geneID, x, y, Dict, debug=False):
    # pdb.set_trace()
    gID = geneID[7:]
    oldPars, oldKids = ({}, {})
    newNames, oldNames = ([], [])
    if x == 0:
        oldPars = {k: Dict[k] for k in Dict.keys() if k[6 : 6 + len(gID)].lower() == gID.lower()}
    elif x == 1:
        oldKids = {k: Dict[k] for k in Dict.keys() if k[-9 : -9 + len(gID)].lower() == gID.lower()}
    elif x == 2:
        oldPars = {k: Dict[k] for k in Dict.keys() if k[6 : 6 + len(gID)].lower() == gID.lower()}
        oldKids = {k: Dict[k] for k in Dict.keys() if k[-9 : -9 + len(gID)].lower() == gID.lower()}
    for par in oldPars:
        Dict[capID(capID(par,6), y)] = Dict.pop(par)
        newNames.append(capID(capID(par,6), y))
        oldNames.append(par)
    for kid in oldKids:
        Dict[capID(capID(kid,-9), y)] = Dict.pop(kid)
        newNames.append(capID(capID(kid,-9), y))
        oldNames.append(kid)
    return(Dict)

def modMeta(splits, mergeds, x):
    def mod(Dict): return({capID(gene,x): Dict[gene] for gene in Dict})
    return(mod(splits), mod(mergeds))

def addGenes(newSplits, oldSplits, newMergeds, oldMergeds, verbose=False):
    oldSplits, newSplits = resolveConflicts(newSplits, oldSplits, verbose)
    oldMergeds, newMergeds = resolveConflicts(newMergeds, oldMergeds, verbose)
    splits = dict(oldSplits, **newSplits)
    mergeds = dict(oldMergeds, **newMergeds)
    return(splits, mergeds)

def resolveConflicts(newGenes, Dict, verbose=False):
    def metaSum(metaStr): return(int(metaStr[1]) + int(metaStr[3]) + int(metaStr[5]))
    conflicts = {k: newGenes[k] for k in newGenes.keys() if k[7:].lower() in [x[7:].lower() for x in Dict.keys()]}
    newChilds = {k: newGenes[k] for k in newGenes.keys() if k[7:].lower() not in [x[7:].lower() for x in Dict.keys()]}
    for conflict in conflicts:
        gID = conflict[6:]
        metaID = conflict[:6]
        conflicted = {k: Dict[k] for k in Dict.keys() if k[6:].lower() == gID.lower()}
        uppers = [k for k in [6, 13] if conflict[k].isupper()]
        for con in conflicted:
            uppers = [k for k in [6, 13] if con[k].isupper()]
            if metaSum(metaID) <  metaSum(con):
                metaID = con[6:]
        newID = metaID + conflict[6:]
        for cap in uppers:
            newID = capID(newID, cap)
        if verbose: print(f"resolved {conflict} to {newID}")
        # pdb.set_trace()
    return(Dict, newChilds)

def update(parGene, childGenes, mergeds, splits, parents):
    if type(parents) != list:
        ref_letter = MetaKey2[parents[7:][0]]
        parGene, childGenes = modMeta(parGene, childGenes, ref_letter)    
        splits = modDict(parents, 2, ref_letter, splits)
        splits, mergeds = addGenes(childGenes, splits, parGene, mergeds, args.v)
        mergeds = modDict(parents, 0, ref_letter, mergeds)
    else:
        for par in parents:
            ref_letter = MetaKey2[par[7:][0]]
            parGene, childGenes = modMeta(parGene, childGenes, ref_letter)
            mergeds = modDict(par, 0, ref_letter, mergeds)
            splits = modDict(par, 2, ref_letter, splits)
        childGenes = {k: childGenes[k] for k in childGenes.keys() if k.lower() not in [x.lower() for x in splits.keys()]}
        splits, mergeds = addGenes(childGenes, splits, parGene, mergeds)
    return(splits, mergeds) 

def writeOutput(split_gene_key, file_path, coord_dict, to_coords):
    with open(file_path, 'w') as outfile:
        for gene in split_gene_key:
            refID = gene[6]
            try:
                outlist = [str(f) for f in sorted([int(x) for x in coord_dict[refID][gene]])]
            except KeyError:
                outlist = [-9, -9, -9]
            outlist.append(gene)
            entry_list = []
            num_entries = len(SGK2[gene])
            for entry in SGK2[gene]:
                if any(k[6] == to_coords for k in entry): # Check if entry has any genes from current ref being converted to 
                    entry_list.append(",".join(entry))
                    entry_str = ",".join(entry[1:])
                    if outlist != [-9, -9, -9, gene]:
                        outfile.write("\t".join(outlist + [entry[0]]) + "\t" + entry_str + "\t" + str(num_entries) + "\n")
                # entry_str = ";".join(entry_list)
            
    return()

def newNames(split_gene_key, coord_dict, to_coords, outfile):
    mergeds = {}
    splits = {}
    parent_list = []
    child_list = []
    debug = False
    cc, dd, ee = (0, 0, 0)
    for parent in split_gene_key:
        # if "b007365" in parent: debug = True
        Child = False
        Grandparent = False
        new_kin_set = []
        for prog_set in split_gene_key[parent]:
            if any(k in parent_list for k in prog_set): # Current child is previously a parent
                Grandparent = True
                cc += 1
                childPars = [k for k in prog_set if k in parent_list]
                new_kin_set.append([capID(x,7) if x in childPars else x for x in prog_set])
            if parent in child_list: #Current parent is previous child
                Child = True
                dd += 1
                print(f"Nested event found for child/parent {[k for k in prog_set if k in parent_list]} and (grand)parent {parent}")
            if any(k in child_list for k in prog_set):
                ee += 1
                print(f"Child gene(s) {[k for k in prog_set if k in child_list]} has multiple parents in addition to {parent}")
        if not Grandparent and not Child:
            parGene, childGenes = newGenes(parent, split_gene_key[parent], coord_dict, to_coords)
            splits, mergeds = addGenes(childGenes, splits, parGene, mergeds)
        if Grandparent: # Something is going on with redundant calls to newGenes or update when gene is grandparent and child
            if debug: pdb.set_trace()
            parGene, childGenes = newGenes(parent, new_kin_set, coord_dict, to_coords)
            splits, mergeds = update(parGene, childGenes, mergeds, splits, childPars)
            if debug: pdb.set_trace()
        if Child:
            if debug: pdb.set_trace()
            parGene, childGenes = newGenes(capID(parent,7), split_gene_key[parent], coord_dict, to_coords)
            splits, mergeds = update(parGene, childGenes, mergeds, splits, parent)
            if debug: pdb.set_trace()
        for prog_set in split_gene_key[parent]:
            child_list = child_list + prog_set
        parent_list.append(parent)

    with open(outfile.replace(".bed", ".merged.bed"), 'w') as outfile1:
        for gene in mergeds:
            coords = "\t".join([str(k) for k in sorted(mergeds[gene])])
            outfile1.write(f"{coords}\t{gene}\n")
    with open(outfile.replace(".bed", ".split.bed"), 'w') as outfile2:
        for gene in splits:
            coords = "\t".join([str(k) for k in sorted(splits[gene])])
            outfile2.write(f"{coords}\t{gene}\n")

    newgenes = []
    with open(outfile.replace(".bed", ".split.bed"), 'r') as file:
        for line in file:
            line = line.strip("\n").split("\t")
            newgenes.append(line[-1][6:-2])
    ff = 0
    for i, gene in enumerate(newgenes):
        switched = gene[7:] + gene[:7]
        for cgene in newgenes[i+1:]:
            if cgene.lower() == switched.lower(): 
                print(cgene, switched)
                ff += 1

    print(cc, dd, ee, ff) 
    return()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "This script takes a series of split/merged candidates and produces a single pseudo-annotation file containing coordinates of each gene along with their corresponding split/merged genes.  The coordinates within the bed files (input in -b) should all be in terms of a common reference as accomplished with convert_coords_wExons.py")
    parser.add_argument('-s', type=str, metavar='split-gene-key', required=True, help='comma separated file where first entry is the parent/merged gene and all subsequent entries are corresponding genes from split model in alternative reference')
    parser.add_argument('-b', type=str, metavar='bed_file_list', required=True, help='comma separated string listing bed files to combine into a pseudo-annotation.  All entries should have coordinates in terms of a common reference')
    parser.add_argument('-c', type=str, metavar='to_coords', required=False, default="-9", help = "Coordinates you are converting to: 1, 4, or 8 for B73, W22, and PH207, respectively")
    parser.add_argument('-o', type=str, metavar='output_file', required=True, help = "output bed file ")
    parser.add_argument('--newNames', action="store_true", help="use new naming scheme for bedfile")
    parser.add_argument('-v', action="store_true", help="verbose")
    args = parser.parse_args()

    bed_files = args.b.split(',')
    coord_dict = makeBedDict(bed_files, args.v)    
    SGK, SGK2 = makeSplitGeneKeyDict(args.s)
    if not args.newNames:
        writeOutput(SGK2, args.o, coord_dict, args.c)
    else:
        if args.c == "-9":
            print("Must specify -c when using --newNames")
        else:
            newNames(SGK, coord_dict, args.c, args.o)
 
    # need to assert that split genes fall within boundaries proposed by merged gene

    