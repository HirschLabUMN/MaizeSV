# Genome STRiP: Preparing a reference sequence

Create an empty folder and create/place the files below into this directory.  Italic text below is taken directly from: http://software.broadinstitute.org/software/genomestrip/.  The metadata bundles contain:

#### reference.fasta
_The reference genome sequence (and accompanying index file in reference.fasta.fai)._

You must also index the fasta using bwa index.  However, you have to use the bwa version that ships with Genome STRiP.

    ~/software/svtoolkit/bwa/bwa index -a bwtsw B73_chr1-10.fasta

#### reference.ploidymap.txt
_This file specifies the gender-dependent normal ploidy of each segment of the reference genome and is used in structural variation calling on the sex chromosomes. In humans, this file flags the Y chromosome as being present only in males and indicates that the X chromosome differs in normal copy number based on the sex of the individual (except for the pseudo-autosomal regions). The structure of this file is current human-centric._

For a hermaphroditic species, these files simply contain the following single line:

    * * * * 2

#### reference.svmask.fasta

_A genome alignability mask (in indexed fasta format) that is used to deterimine whether reference sequence should be considered to be uniquely alignable.
This mask is computed based on the reference genome and a default k-mer length and indicates whether a k-mer centered on each position of the reference genome has a unique mapping to the reference. The k-mer length is defaulted in each reference metadata bundle, but this mask can be overridden if necessary for specific data sets._

_This file is usually present in each reference metadata bundle as a symbolic link to a file whose name specifies the k-mer length used to compute the mask._

Genome STRiP has an utility (ComputeGenomeMask) that will create the svmask.fasta file, which is implemented in *GSTRiP_ComputeGenomeMask.sh*.  This script parallelizes the mask computation across chromosomes, concatenates results, and indexes the concatenated fasta file using Samtools.  The full path(s) for the reference(s) that you want to compute the svmask.fasta for need to be provided in a file (one path per line).  The path to this file should be specified as REF_PATH within the shell script.  There are numerous paths within this script that will need to be modified for the user.

#### reference.gcmask.fasta

_A genome mask (in indexed fasta format) that is used for gc-bias estimation in each sequencing library. This mask excludes region of the genome that have variable ploidy (e.g. sex chromosomes), regions that are in annotated segmental duplications, regions that have repeats flagged by repeat masker, regions annotated as containing polymorphic CNVs in DGV and miscellaneous reference contigs. In older versions of Genome STRiP, this file was called the "cn2mask" because it selects regions of the genome that are likely to be two copies in most people.
There are also an accompanying index file (reference.gcmask.fasta.fai), a file with statistics about the number of masked and unmasked bases (reference.gcmask.fasta.stats), and for convenience a bed file with the equivalent mask information (reference.gcmask.fasta.bed)._

The GC mask is created via another Genome STRiP utility (SVPreprocess).  This is a relatively long and error-prone process that may need to be resubmitted multiple times.  The SVPreprocess.q script that ships with Genome STRiP fails to find the correct path for an additional script necessary to estimate gender from the sequence data.  Although we don’t need to estimate gender, there is no way to prevent program from trying.  A corrected version of this script is SVPreprocess1.q in _scripts_, but you can also change for yourself by changing:

    val plotReadDepthScriptPath = "/home/hirschc1/pmonnaha/software/svtoolkit/R/metadata/plot_chr_vs_chr_readdepth.R"
within the SVPreprocess.q script to point to your copy of this R script.  Save this script as SVPreprocess1.q alongside the original.

#### Other files

There are a number of other files in the tutorial linked above, but these are not necessary for running the program.  They deal primarily with issues associated with sex chromosomes.





