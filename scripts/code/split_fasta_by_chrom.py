import sys
from Bio import SeqIO

out_dir = "/" + sys.argv[2].strip("/") + "/"
with open(sys.argv[1]) as original_fasta:
    records = SeqIO.parse(original_fasta, 'fasta')
    for record in records:
    	print out_dir + record.id +".fa"
        SeqIO.write(record, out_dir + str(record.id) + ".fa", 'fasta')