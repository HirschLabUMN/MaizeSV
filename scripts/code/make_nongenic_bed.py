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
args = parser.parse_args()


first_gene = True
old_stop = 0

# Begin looping over lines in the gff file
with open(args.gff, 'r') as gff:
	for i, line in enumerate(gff):
		if line[0] != "#": # Ignore all lines in the header of the gff
			if line.split()[1] != "wareLab": # This is a catch for the B73 gff because it has a wierd first line
				line = line.strip("\n").split()
				chrom = line[0]
				start = line[3] # Lower boundary of entry
				stop = line[4] # Upper boundary of entry

				# Initialize current_chrom to catch when we move on to a different chromosome
				if first_gene is True:
					current_chrom = chrom 
				
				# this catches a gene entry on same chomosome as before
				if line[2] == "gene" and chrom == current_chrom: 
					if int(start) - int(old_stop) >= args.b * 2 : # If buffer regions overlap, then we consider the space between genes to be GENIC.
						if first_gene is True:
							print(chrom + "\t0\t" + str(int(start) - args.b) + "\n")
							first_gene = False
						else:
							print(chrom + "\t" + str(int(old_stop) + args.b) + "\t" + str(int(start) - args.b) + "\n")
				
				# This catches a gene entry after moving on to new chromosome
				elif line[2] == "gene" and chrom != current_chrom: 
					print(current_chrom + "\t" + str(int(old_stop) + args.b) + "\t999999999\n") # the final non-genic region is the end of the last gene stop point all the way to the end of the chromosome (999999999)
					if int(start) - args.b > 0: # Does buffer region of first gene on new chromosome extend off the beginning of the chromosome?
						print(chrom + "\t0\t" + str(int(start) - args.b) + "\n")
					else:
						print(chrom + "\t0\t" + str(start) + "\n")
					current_chrom = chrom
				old_stop = stop # Store upper boundary of current gene 




