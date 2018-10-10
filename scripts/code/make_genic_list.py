#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: This script creates a bed file of NON-genic regions to be excluded in structural variant discovery
 Takes the following arguments:
    -gff :  Full path to gff file containing gene locations (can contain other elements as well...these will just be ignored)
	-b   :  Buffer on either side of gene boundary.  E.g Gene1 ends at 100bp and Gene2 starts at 500bp.  If b=10, then the non-genic region will be 110-490')
"""

# Import necessary arguments
import argparse

# Specify arguments to be read from the command line
parser = argparse.ArgumentParser()
parser.add_argument('-gff', type=str, metavar='gff_path', required=True, help='path to gff file')
parser.add_argument('-b', type=int, metavar='buffer', default=2000, help='Buffer on either side of gene boundary.  E.g Gene1 ends at 100bp and Gene2 starts at 500bp.  If b=10, then the non-genic region will be 110-490')
parser.add_argument('-r', action='store_true', help='restrictTo_mainChroms')
args = parser.parse_args()


def merge_intervals(intervals):
	# taken from https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
	sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
	merged = []

	for higher in sorted_by_lower_bound:
	    if not merged:
	        merged.append(higher)
	    else:
	        lower = merged[-1]
	        # test for intersection between lower and higher:
	        # we know via sorting that lower[0] <= higher[0]
	        if higher[0] <= lower[1]:
	            upper_bound = max(lower[1], higher[1])
	            merged[-1] = (lower[0], upper_bound)  # replace by merged interval
	        else:
	            merged.append(higher)
	return merged

old_chrom = "99"
included_chrs = ['1','2','3','4','5','6','7','8','9','10','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr01','chr02','chr03','chr04','chr05','chr06','chr07','chr08','chr09','chr10']
Unmerged_Intervals = {}
# Begin looping over lines in the gff file
with open(args.gff, 'r') as gff:
	for i, line in enumerate(gff):
		if line[0] != "#": # Ignore all lines in the header of the gff
			if line.split()[1] != "wareLab": # This is a catch for the B73 gff because it has a wierd first line
				line = line.strip("\n").split()
				chrom = line[0].replace("M","chr")
				start = int(line[3]) # Lower boundary of entry
				stop = int(line[4]) # Upper boundary of entry				
				if line[2] == "gene":
					if start - args.b < 0: 
							start = 0
					else:
						start = start - args.b
					if chrom != old_chrom:
						if args.r is True and chrom in included_chrs:
							Unmerged_Intervals[chrom] = [(start,stop)]
						elif args.r is False:
							Unmerged_Intervals[chrom] = [(start,stop)]
						old_chrom = chrom
					else:
						if args.r is True and chrom in included_chrs:
							Unmerged_Intervals[chrom].append((start,stop))
						elif args.r is False:
							Unmerged_Intervals[chrom].append((start,stop))


for chrom in Unmerged_Intervals:
	Merged_Intervals = merge_intervals(Unmerged_Intervals[chrom])
	for i in Merged_Intervals:
		print(f"{chrom}:{i[0]}-{i[1]}")

