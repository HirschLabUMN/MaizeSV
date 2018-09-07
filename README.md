# Structural Variation Pipeline

#### Summary

Structural variant discovery using whole genome short-read sequence data.  Reads are preprocessed with Sickle and CutAdapt, then mapped to multiple refernce genomes using SpeedSeq. Currently, structural variants are identified with LUMPY and Genome STRiP. Custom scripts to generate consensus call set across reference genomes and softwares.  Shell scripts within the ./jobs directory are designed to implement the softwares and scripts on an HPC using an SGE job submission system. 

## Requirements

### Python
Python 2.7 is required for SVtools (a component of the Lumpy pipeline)

To create a virtual environment with anaconda using python 2.7, do:

    conda create -n py27 python=2.7 anaconda

Activate this environment prior to using any SVtools or other py27-dependent softwares with:

    source activate py27

You will also need the following python modules, which can be installed with pip or conda:

    pip install biopython
    pip install --user --upgrade cutadapt
    conda install -c bioconda svtools
  
### _Speedseq_ <https://github.com/hall-lab/speedseq>
Also installs _sambamba_ and _samblaster_.

    module load cmake
    git clone --recursive <https://github.com/hall-lab/speedseq>
    cd speedseq
    make

If the installation fails, try installing just the necessary components of speedseq, which is align and sv.  Installation of speedseq is modular, so these components can be installed with individual calls to make.
E.g. 

    make align
from within the speedseq directory.

For installing _speedseq sv_, first do:
    
    module load root
    source /panfs/roc/msisoft/root/6.06.06/bin/thisroot.sh

### _RepeatMasker_

Both Genome STRiP and Lumpy (svtools classify) require a list of mobile element insertion sites, which are most easily produced via RepeatMasker <http://www.repeatmasker.org/>

Implemented with ./scripts/jobs/RepeatMasker.sh.  

### _Lumpy_ <https://github.com/arq5x/lumpy-sv>

We implement Lumpy on a per-individual basis via the _speedseq_ 'sv' command

For installing speedseq sv, do:

    module load root
    source /panfs/roc/msisoft/root/5.34.32/bin/thisroot.sh
    make sv
    make cnvnator

NOTE:  It is essential that you use this specific version of ROOT for both compiling 

Once the installations are complete, you must modify the file ./speedseq/bin/speedseq.config in a couple of ways.  First add the line:

    source /panfs/roc/msisoft/root/5.34.32/bin/thisroot.sh

under the #CNVnator heading.

Next, add the virtual environment python path to to where it says “PYTHON=“.  For me, the python path looks like:

    home/hirschc1/pmonnaha/anaconda3/envs/py27/bin/python2.7
    

## Fasta preparation

Several subsequent steps in the pipeline work more smoothly if all extraneous scaffolds are removed from the reference sequences.  Use **filter_fasta.py** to achieve this:

To install all necessary components to run **filter_fasta.py**, perform the following steps (taken from <https://www.biostars.org/p/157811/>)

1.Ensure you have python and biopython installed. Type in your terminal:  

    source activate py27
    python -c "Import Bio"
    echo $?
    sudo pip install biopython

2.Run the program as follows: 

    chmod +x filter_fasta.py
    python filter_fasta.py your_ids.txt IN.fasta OUT.fasta

where: 

  * your_ids.txt --> A file containing the identifiers including ">" you want to EXCLUDE, one identifier per line; like: 

``>scaffold_999`` 

  * IN.fasta --> your original fasta file

  * OUT.fasta --> your output file
<br />

Files ending with **extraContigs** in _/scripts/accessorry/_ contain the extra contig ID’s for references from the following paths within _/home/maize/_:

    B73  : ./shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa
    PH207: ./shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa
    W22  : ./shared/databases/genomes/Zea_mays/W22/W22__Ver12.genome.normalized.fasta   
    PHB47: ./sna/PHB47/Zea_mays_var_PHB47.mainGenome.fasta  

After removing extraneous scaffolds, reference fastas need to be indexed using _bwa_:

    bwa index <reference_fasta>



## Fastq Pre-processing

### Quality Assessment — _Fastqc_

To generate _fastqc_ commands:

    find <fastqdir> -name "*fastq*" -o -name "*fq*" -print | xargs -I {} \
    echo "fastqc --noextract -t number_of_threads -o outdir" {} > command_file
  
You can use **fastqc.sh** in _/scripts/jobs/_ to run these commands in parallel via a task array.  If you are unfamiliar with the different options for parallelizing job submission, see the discussion below in **Notes**.  You will also need to open **fastqc.sh** and modify paths where indicated.
<br />

### Adapter and Quality Trimming — _CutAdapt_ and _Sickle_

We use _cutadapt_ and _sickle_ to trim adapters and low quality bases, respectively, from reads.  Both of these softwares need to be downloaded and installed.  Although _cutadapt_ is available as load-able module on MSI, the correct version is not available.

Install _cutadapt_ using pip

    pip install --user --upgrade cutadapt

Install _sickle_ by cloning the git repository and typing ‘make’ from within the directory

    git clone <https://github.com/najoshi/sickle.git>
    cd sickle
    make

To generate commands for cutadapt and sickle, use the script _./scripts/code/_**Generate_CutAdapt_commands.py**.  The use of this script differs depending on whether the paired-end fastq data is in an interleaved format (both forward and reverse reads are interleaved in a single file) versus in separate files.  For interleaved data, simply do:

    python Generate_CutAdapt_commands.py -c path_to_cutadapt \
    -s path_to_sickle \
    -r read_path_file \
    -o output_directory >> commands_file

whereas for data in separate files do:

    python Generate_CutAdapt_commands.py -c path_to_cutadapt \
    -s path_to_sickle \
    -r read1_path_file \
    -r2 read2_path_file \
    -o output_directory >> commands_file

where the _read1_path_file_ contains the forward read files and _read2_path_file_ contains the paths to the reverse reads.  

**NOTE**: the ‘>>’ in the above commands will append the output of each run to the commands_file, which is useful for combining commands for different sets of interleaved/separate fastq files.  However, be sure the commands_file doesn’t already exist, or you may potentially append on top of a bunch of old commands that you don’t want to include.

Run commands as a task array with script _./scripts/jobs/_**Cutadapt.sh**.  Open **Cutadapt.sh** and modify the paths where indicated prior to running.


## Mapping Fastq Files

#### Alignment - _Speedseq_

We use _speedseq_ for mapping, position-sorting, and extraction of split and discordant paired reads.  Under the hood, _speedseq_ more efficiently parallelizes _bwa mem_ and pipes the output directly to _sambamba_ and _samblaster_.  _Sambamba_ is a sam/bam manipulation software with analogous functionality to _samtools_, like sorting, merging, etc., yet is much faster.  _Samblaster_ efficiently marks duplicate reads and simultaneously extracts discordant and split reads.

Install _speedseq_ (also installs _sambamba_ and _samblaster,_ automatically) with:

    module load cmake
    git clone --recursive <https://github.com/hall-lab/speedseq>
    cd speedseq
    make

If the installation fails, try installing just the necessary components of speedseq, which is align and sv.  Installation of speedseq is modular, so these components can be installed with individual calls to make.
E.g. 

    make align
from within the speedseq directory.

For installing _speedseq sv_, first do:
    
    module load root
    source /panfs/roc/msisoft/root/6.06.06/bin/thisroot.sh

Generate _speedseq_ commands with _./scripts/code/_**Generate_SpeedSeq_commands.py**

    python Generate_SpeedSeq_commands.py -f fastq_directory \
    -s Sample_fastq_key  \
    -o output_directory \
    -r reference_path_file \
    -c number_of_cores \
    -m number_of_Gb_memory > speedseq_command_file

  
Given the large number of nodes needed to complete all mapping commands, it is most efficient to use either the ‘large’ or ‘widest’ queues on mesabi.  To use the ‘widest’ queue, you must request special permission by emailing help@msi.umn.edu.  

For the ‘large’ queue, the maximum number of nodes is 48 and the maximum time limit is 24 hrs.  Trials with _speedseq_  suggest ~8hrs is a reasonable expectation for a single mapping job.  Therefore, ~3 jobs per node ought to complete within the time limit, resulting in (48 nodes * 3 jobs/node) 144 total jobs per run in the large queue.  We can split the _speedseq_command_file_ into multiple files each containing 144 lines using: 

    split -l 144 speedseq_command_file -a 1 -d speedseq_commands_

This will produce a number of files named _speedseq_commands_X_ where X will be replaced by a numerical index.  The job script ./scripts/jobs/ **Speedseq_large.sh,** implemented as a task array, can be used to submit a job to the ‘large’ queue for each of the subsetted command files.  Open **Speedseq_large.sh** and modify the paths where indicated.
<br />

#### Merging BAMs — _sambamba_

The script  **Generate_MergeBAMs_commands.py**  in _./scripts/code/_ can be used to generate commands to merge bams with _sambamba:_

    python Generate_MergeBAMs_commands.py -k _sample_fastq_key_  \
    -o _output_directory_ \
    -r reference_path_file \
    -c _number_of_cores_  > _sambamba_command_file_

Run _sambamba_  commands as task array implemented in ./scripts/jobs/ **Sambamba_MergeBAMs.sh**.  Be sure to open **Sambamba_MergeBAMs.sh** __ and modify paths where necessary.  ALSO, be sure to modify the ‘ppn’ field to reflect the number of cores specified with **Generate_MergeBAMs_commands.py** and also adjust the ‘mem’ field appropriately. 

  
## SV Discovery

The end goal for each of the softwares is to produce a multi-sample VCF for each of the reference assemblies that were mapped to.  We restrict each program to perform SV discovery within gene regions to avoid the massive number of spurious calls that would result from the high TE content outside of genic regions.

To create the bed files that will contain the non-genic regions to be excluded for each reference, use the script ./scripts/code/make_nongenic_bed.py:

    python make_nongenic_bed.py -gff gff_file -b buffer_size > {REF}.NonGenicMask.bed

The -b option specifies a buffer on either side of gene boundary.  E.g Gene1 ends at 100bp and Gene2 starts at 500bp.  If b=10, then the non-genic region will be 110-490.  Current analysis are being done with -b 2000.

### Lumpy

The team behind Lumpy has created several software packages (speedseq, svtyper, and svtools) all of which are an encouraged or necessary part of a pipeline to go from read data to a final multi-sample VCF.

We implement Lumpy on a per-individual basis via the speedseq sv command, which will subsequently run svtyper and CNVnator, including all results in a single output file per individual.  Subsequent steps (implemented ) latter steps are important for merging variant calls across individuals using svtools.

NOTE: See instructions above (Requirements) for installing ROOT prior to proceeding

Prior to running speedseq sv, you must also split each fasta by chromosome and label them exactly as the sequence name appears in the original fasta (filenames must end with ".fa").  For each reference, create a directory containing ONLY these chromosome fastas.  Make one copy of the ./speedseq/bin/speedseq.config file per reference genome and provide the directory containing the chromosome fastas after CNVNATOR_CHROMS_DIR=.  Name these config files speedseq.REF.config where the REF name matches that used in the reference_path_file.

Note: Make certain to run speedseq sv using the following options: -v -d -P -g -k option as subsequent steps will utilize CNVnator files in the temporary directories, assume that SVTyper has been run and require LUMPY's probability curves.

To generate commands for speedseq sv, use the script ./scripts/code/Generate_SpeedSeqSV_commands.py:

    python Generate_SpeedSeqSV_commands.py -b bam_directory \
    -r reference_path_file \
    -o output_directory \
    -c number_of_cores \
    -s speedseq_directory

The reference_path_file is the same as before, but with an additional column containing the path to the non-genic bed files created with make_nongenic_bed.py.

Once all files have run through speedseq sv, we can start the process of merging variants across individuals to create the final multi-sample VCF.  These steps will be performed using a number of functions within svtools.  svtools is single-threaded so parallelization is implemented via GNU parallel in all corresponding shell scripts.

You can download and install svtools plus all requirements with:

    source activate py27
    conda install -c bioconda svtools

NOTE: You must activate this virtual environment every time you wish to run svtools.

After everything installs successfully, the first step is to use svtools lsort, which will combine all of the individual VCF files produced from speedseq sv into a single, sorted VCF followed by svtools lmerge, which will merge overlapping SV calls.

These two processes are implemented in ./scripts/jobs/SVtools_SortAndMerge.sh as follows:

    qsub SVtools_SortAndMerge.sh -F "input_vcf_dir output_vcf temp_dir batch_size percent_slop” \
    -l mem=20gb,walltime=24:00:00 \
    -q mesabi \
    -o stdout_file \
    -e stderr_file

All arguments are required, must follow correct order, and must be space-delimited within quotations following the -F flag.  You can change the requested memory and walltime using the -l flag.  

_input_vcf_dir_ is directory containing per-individual VCFs.  If you have called SVs against multiple reference genomes, create separate directories for each reference and call the script once per input directory.  output_vcf is the full path to output filed cannot already exist.  temp_dir is where to place temporary files during the sorting step.  batch_size is the number of vcfs to sort at a time, which determines amount of memory required.  percent_slop is the percent increase of the breakpoint confidence interval both up and down stream, which determines the aggressiveness of the merging process. 

The next step is to run the svtools genotype and svtools copynumber.  Prior to this, you must create the coordinates_file needed for CNVnator for each of the merged vcf files.

    zcat merged_ref.vcf.gz > merged_ref.vcf
    source activate py27
    create_coordinates -i merged_ref.vcf -o ref_coordinates_file
_create_coordinates_ is a script that ships with svtools and should install as an executable in your /usr/bin/ directory.  It will be available after activating the python 2.7 virtual environment in which svtools was installed.

To generate commands for svtools genotype and svtools copynumber use the script in ./scripts/code/Generate_SVtools_Genotype_commands.py

    python Generate_SVtools_Genotype_commands.py \
    -b bam_directory_used_in_speedseq_sv \
    -o output_directory \
    -v merged_ref.vcf.gz \
    -c ref_coordinates_file \
    -s path_to_speedseq_software_directory \
    -w window_size >> svtools_genotype_commands_file
You will need to run this once for each merged_ref.vcf.  These commands can be run using ./scripts/jobs/SVtools_Genotype.sh, which uses GNU parallel.  They can also be split and run as multiple GNU parallel jobs via a task array as we did with Speedseq_large.sh.

NOTE: The ‘>>’ will append commands to a file, so if an older file already exists, you should delete the older file before running this script.

The final step in the Lumpy/svtools pipeline involves pasting the individual genotyped vcfs that include information from CNVnator, further prune the combined VCF for redundant SV calls, and classify remaining calls as valid or artificial.  Prune takes a clustering approach based on a specified evaluation parameter, and prunes sites in the same cluster that are within a given distance.

Use the script ./scripts/jobs/SVtools_PastePruneClassify.sh:

    qsub SVtools_PastePruneClassify.sh -F "vcf_dir merged.vcf out_prefix dist eval_param te.bed”
    -l mem=20gb,walltime=24:00:00
    -o stdout_file \
    -e stderr_file

vcf_dir is the input directory containing per-individual VCFs.  If you have called SVs against multiple reference genomes, create separate directories for each reference and call the script once per input directory.  merged.vcf is the merged vcf that was used in svtools genotype/copynumber.  out_prefix is the prefix to use for output vcf.  These files cannot already exist.  dist is the max separation distance (bp) of adjacent loci in a cluster [50].  eval_param is the evaluating parameter for choosing best bedpe in a cluster(e.g. af=AlleleFrequency default:af).  Currently, it seems that “af” is the only option for eval_param.  te.bed is a gzipped bed file from repeat masker or elsewhere that specifies te locations.

### Genome STRiP

Prior to running Genome STRiP, you have to create a “MetaData Bundle” for each reference.  Create an empty folder and create/place the files below into this directory.  Italic text below is taken directly from: http://software.broadinstitute.org/software/genomestrip/.  The metadata bundles contain:

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

Wrote job script to run Genome STRiP utility (ComputeGenomeMask).  Parallelizes by scaffold.  Change file to point to file containing paths to references that you want to perform operation on.  One path per line.

### Notes

* * *

#### Parallelization

There are two options for parallelizing jobs on the MSI cluster:

  1. GNU parallel 

  2. Task Arrays

MSI has good help documentation for each approach.  For task arrays, see <https://www.msi.umn.edu/support/faq/how-do-i-use-job-array>.  For GNU parallel, see <https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel>

The primary difference is that GNU parallel, for the most part, can only handle single-threaded jobs, whereas task arrays can be used for multi-threaded jobs.  Task arrays also allow for simpler resubmission of failed jobs.

Compare the two scripts **fastqc_GNUparallel.sh** and **fastqc_TaskArray.sh** in _./scripts/jobs/_ to see how job specification differs between the two approaches.  To submit these jobs, you would simply use:

    qsub fastqc_GNUparallel.sh

for the former, and:

    qsub -t 0-500 fastqc_TaskArray.sh

for the latter, replacing ‘500’ for the total number of commands in your command file.  Failed jobs can be resubmitted with ``qsub -t job_number``

  

