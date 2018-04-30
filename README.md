# WiDiv Structural Variation

Summary

* * *

Structural variant discovery using whole genome short-read sequence data on 500 inbred lines from the Wisconsin Diversity Panel (WiDiv).  For robust structural variant discovery, we take the consensus variant calls across:

  


  1. Three state-of-the-art SV discovery tools

    * Lumpy 

    * MetaSV

    * Genome STRiP

  2. At least 4 high-quality reference genomes

    * B73, W22, PH207, and PHB47




  


Sequence data comes from 3 sources:

  1. JGI: An initial subset of ~60 lines

  2. Mark Mikel (MM): High depth sequence data corresponding to 8 lines

    * B73, G84, G35, J40, G39, LH82, G47, and PH207

  3. Novagene: All remaining samples




  


Fasta preparation

* * *

Several subsequent steps in the pipeline work more smoothly if all extraneous scaffolds are removed from the reference sequences.

  


Use filter_fasta.py to achieve this:

python filter_fasta.py Ref_extraContigs.txt IN.fasta OUT.fasta

  


To install all necessary components to run  **filter_fasta.py,** perform the following steps (taken from <https://www.biostars.org/p/157811/>)

  


1) Ensure you have python and biopython installed. Type in your terminal: 

python -c "import Bio"

  


echo $?

If "0" displayed, pass to step 2. In case error message appears, install it by easy_install or pip: 

sudo easy_install -f <http://biopython.org/DIST/> biopython

or 

sudo pip install biopython

In case of not having pip installed, type: 

sudo apt-get install pip

  


3) Run the program as follows: 

chmod +x filter_fasta.py

  


python filter_fasta.py your_ids.txt IN.fasta OUT.fasta

where: 

  * your_ids.txt --> A file containing the identifiers including ">" you want to EXCLUDE, one identifier per line; like: 




>hhh 

>ghag 

>MAMND

  * IN.fasta --> your original fasta file

  * OUT.fasta --> your output file




  


Files ending with **extraContigs** in /scripts/accessorry/ contain the extra contig ID’s for references from the following paths within /home/maize/:

B73  : ./shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa

PH207: ./shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa

W22  : ./shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta   

PHB47: ./sna/PHB47/Zea_mays_var_PHB47.mainGenome.fasta

  


After removing extraneous scaffolds, reference fastas need to be indexed using bwa:

bwa index <reference_fasta>

  


Fastq Pre-processing

* * *

Assess quality of fastq data via fastqc as implemented in .  

  


The sequence data from MM (8 samples listed below) need to be handled slightly differently

  * B73, G84, G35, J40, G39, LH82, G47, and PH207




  


These samples need to be handled slightly differently.  Concatenated all fastq files following adapter and quality trimming
