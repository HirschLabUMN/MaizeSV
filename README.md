# Structural Variation Pipeline

#### Summary

Structural variant discovery using whole genome short-read sequence data.  Reads are preprocessed with Sickle and CutAdapt, then mapped to multiple refernce genomes using SpeedSeq. Currently, structural variants are identified with LUMPY and Genome STRiP. Custom scripts to generate consensus call set across reference genomes and softwares.  Shell scripts within the ./jobs directory are designed to implement the softwares and scripts on an HPC using an SGE job submission system. 

#### Important Notes

The "< >" notation indicates a parameter, path, etc. that needs to be specified by the user when running the code/command.  For instance, the user would specify the path to the reference fasta file in the following command WITHOUT including the "< >".

        bwa index <reference_fasta>

All of the scripts in the ./scripts/code/ directory use python 3.6.  Load this on the cluster by:

    module load python/3.6.3

At the bottom of the README, you can find a list of helpful commands at various stages in the pipeline.

## Requirements

### Python

Python 2.7 is required for SVtools (a component of the Lumpy pipeline)

To create a virtual environment with anaconda using python 2.7, do:

    module load anaconda/1.7.0 
    conda create -n py27 python=2.7 anaconda

Activate this environment prior to using any _SVtools_ or other py27-dependent softwares with:

    source activate py27

You will also need the following python modules, which can be installed with pip or conda:

    pip install biopython
    pip install pysam
    pip install --user --upgrade cutadapt
    pip install scikit-allel
    conda install -c bioconda svtools


### FASTQ Pre-processing
We use _cutadapt_ and _sickle_ to trim adapters and low quality bases, respectively, from reads.  Both of these softwares need to be downloaded and installed.  Although _cutadapt_ is available as load-able module on MSI, the correct version is not available.

Install _cutadapt_ using pip

    pip install --user --upgrade cutadapt

Install _sickle_ by cloning the git repository and typing ‘make’ from within the directory

    git clone https://github.com/najoshi/sickle.git
    cd sickle
    make

  
### _Speedseq_ <https://github.com/hall-lab/speedseq>
Speedseq is used to map fastq files to each reference genome (via the align submodule) and is also used to call SVs with Lumpy (via the sv submodule).  This process will also install _sambamba_ and _samblaster_, which are necessary for mapping and merging of sequence data.

    module load cmake
    module load gcc
    source /panfs/roc/msisoft/root/5.34.32/bin/thisroot.sh
    git clone --recursive https://github.com/hall-lab/speedseq
    cd speedseq

Only the align, sv, and cnvnator modules are necessary for mapping and calling SVs, and these can be compiled individually via the following commands.

    make align
    make sv
    make cnvnator

from within the speedseq directory.

Once the installations are complete, you must modify the file ./speedseq/bin/speedseq.config in a couple of ways.  First add the line:

    source /panfs/roc/msisoft/root/5.34.32/bin/thisroot.sh

under the #CNVnator heading.  NOTE:  It is essential that you use this specific version of ROOT for both compiling and subsequent running of speedseq sv

Next, check that the virtual environment python path that you set up earlier is found where it says “PYTHON=“ in the speedseq.config file.  For me, the python path looks like:

    home/hirschc1/pmonnaha/anaconda3/envs/py27/bin/python2.7


### _RepeatMasker_

    module load repeatmasker

Both Genome STRiP and Lumpy (svtools classify) require a list of mobile element insertion sites, which are most easily produced via RepeatMasker <http://www.repeatmasker.org/>

Implemented with ./scripts/jobs/RepeatMasker.sh.  

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
    python filter_fasta.py <your_ids.txt> <IN.fasta> <OUT.fasta>

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
    echo "fastqc --noextract -t <number_of_threads> -o <outdir>" {} > <command_file>
  
You can use **fastqc.sh** in _/scripts/jobs/_ to run these commands in parallel via a task array.  If you are unfamiliar with the different options for parallelizing job submission, see the discussion below in **Notes**.  You will also need to open **fastqc.sh** and modify paths where indicated.
<br />

### Adapter and Quality Trimming — _CutAdapt_ and _Sickle_

First load python 3.6:

    module load python/3.6.3

To generate commands for cutadapt and sickle, use the script _./scripts/code/_**Generate_CutAdapt_commands.py**.  The use of this script differs depending on whether the paired-end fastq data is in an interleaved format (both forward and reverse reads are interleaved in a single file) versus in separate files.  For interleaved data, simply do:

    python Generate_CutAdapt_commands.py -c <path_to_cutadapt> \
    -s <path_to_sickle> \
    -r <read_path_file> \
    -o <output_directory> >> <commands_file>

whereas for data in separate files do:

    python Generate_CutAdapt_commands.py -c <path_to_cutadapt> \
    -s <path_to_sickle> \
    -r <read1_path_file> \
    -r2 <read2_path_file> \
    -o <output_directory> >> <commands_file>

where the _read1_path_file_ contains the forward read files and _read2_path_file_ contains the paths to the reverse reads.  

**NOTE**: the ‘>>’ in the above commands will append the output of each run to the commands_file, which is useful for combining commands for different sets of interleaved/separate fastq files.  However, be sure the commands_file doesn’t already exist, or you may potentially append on top of a bunch of old commands that you don’t want to include.

Run commands as a task array with script _./scripts/jobs/_**Cutadapt.sh**, via:

    qsub -t <XX>-<YY> Cutadapt.sh -F "<path_to_command_file>"


## Mapping Fastq Files

Before proceeding, I will discuss two accessory files that are of critical importance in several subsequent steps: the *sample_fastq_key* and the *reference_path_file*.  These can be found in ./scripts/accessory/sample_fastq_key.txt and ./scripts/accessory/reference_paths.txt, respectively.  The first file contains the sample names in the first column, and then the fastq file prefixes associated with the sample names in (tab-delimited) subsequent columns.  For example:

    A322    GH17_1001_DDPL02023_HT2FJCCXY_L1    GH17_1001_DDPL02023_HTLWWCCXY_L1    r1001_FDPL190729834-1a_HTWF3DSXX_L1 

indicates that the A322 has 3 sets of fastq files (one file for forward reads and one file for reverse reads per set).  These fastq prefixes should be directly derived from the fastq filenames that are in the raw fastq directory, such that appending '\_1.fq.gz' and '\_2.fq.gz' to the prefix will specify the forward and reverse read filenames, respectively.  Any time that new fastq files are received, this file will need to be updated

The *reference_path_file* is much simpler and will only need to be updated once for a new user.  It will also need to be updated if a new, additional reference genome is used for mapping, but this is not likely to occur.  There are 3 tab-delimited columns in this file.  The first column contains the reference genome name that is used for labelling throughout the pipeline (e.g. B73v4, PH207, PHB47, Mo17, and W22v12).  The second column contains the path to the reference fasta files (including only chromosomes 1-10).  The third column contains the path to the bed files that contain non-genic regions to be excluded.  Creating these bed files is explained below under *SV Discovery*; however, it will likely not be necessary to remake these files.  You should be able to use the NonGenic.bed files in ./scripts/accessory/, as long as you change the paths in the *reference_path_file* to specify where these bed files are located for you.

#### Alignment - _Speedseq_

We use _speedseq_ for mapping, position-sorting, and extraction of split and discordant paired reads.  Under the hood, _speedseq_ more efficiently parallelizes _bwa mem_ and pipes the output directly to _sambamba_ and _samblaster_.  _Sambamba_ is a sam/bam manipulation software with analogous functionality to _samtools_, like sorting, merging, etc., yet is much faster.  _Samblaster_ efficiently marks duplicate reads and simultaneously extracts discordant and split reads.

Generate _speedseq_ commands with _./scripts/code/_**Generate_SpeedSeq_commands.py**

    python Generate_SpeedSeqAlign_commands.py -f <fastq_directory> \
    -k <sample_fastq_key>  \
    -o <output_directory> \
    -r <reference_path_file> \
    -c <number_of_cores> \
    -m <number_of_Gb_memory> \
    -s <path_to_speedseq_directory> > <speedseq_command_file>

From extensive experimentation, I have found that using '-c 18 -m 30', is a good sweet-spot to prevent jobs from exceeding the max memory on a node (and thus failing).  

Given the large number of nodes needed to complete all mapping commands, it is most efficient to use either the ‘large’ or ‘widest’ queues on mesabi.  To use the ‘widest’ queue, you must request special permission by emailing help@msi.umn.edu.  

For the ‘large’ queue, the maximum number of nodes is 48 and the maximum time limit is 24 hrs.  Trials with _speedseq_  suggest ~8hrs is a reasonable expectation for a single mapping job.  Therefore, ~3 jobs per node ought to complete within the time limit, resulting in (48 nodes * 3 jobs/node) 144 total jobs per run in the large queue.  We can split the _speedseq_command_file_ into multiple files each containing 144 lines using: 

    split -l 144 <speedseq_command_file> -a 1 -d speedseq_commands_

This will produce a number of files named _speedseq_commands_X_ where X will be replaced by a numerical index.  The job script ./scripts/jobs/ **Speedseq_large.sh,** implemented as a task array, can be used to submit a job to the ‘large’ queue for each of the subsetted command files.  Open **Speedseq_large.sh** and modify the paths where indicated.  Alternatively, jobs in the command file can be individually submitted as a task array using ./scripts/jobs/Speedseq.sh via

    qsub -t <XX>-<YY> Speedseq.sh -F "<command_file>"

where XX and YY specifies the numeric range (line numbers) of the command file that you wish to submit and *command_file* is the file containing the previously generated commands.
<br />

#### Merging BAMs — _sambamba_

The script  **Generate_MergeBAMs_commands.py**  in _./scripts/code/_ can be used to generate commands to merge bams with _sambamba:_

    python Generate_MergeBAMs_commands.py -k <sample_fastq_key>  \
      -b <bam_directory>
      -o <output_directory> \
      -r <reference_path_file> \
      -c <number_of_cores>  > <sambamba_command_file>

Run _sambamba_  commands as task array implemented in ./scripts/jobs/ **Sambamba_MergeBAMs.sh**.  Be sure to open **Sambamba_MergeBAMs.sh** __ and modify paths where necessary.  ALSO, be sure to modify the ‘ppn’ field to reflect the number of cores specified with **Generate_MergeBAMs_commands.py** and also adjust the ‘mem’ field appropriately. 

You can use the shell script *./scripts/jobs/Sambamba_MergeBAMs.sh* to submit the jobs in the *command_file* as a task array, via:

    qsub -t <XX>-<YY> Sambamba_MergeBAMs.sh -F "<path_to_command_file>"

, where XX and YY correspond to the range of line numbers in the *command_file* that user has specified.

  
## SV Discovery

The end goal for each of the softwares is to produce a multi-sample VCF for each of the reference assemblies that were mapped to.  We restrict each program to perform SV discovery within gene regions to avoid the massive number of spurious calls that would result from the high TE content outside of genic regions.

NOTE: The following bed files have already been created for a buffer size of 2kb, which can be found in ./scripts/accessories.  You only need to run the following scripts if you desire to use a different gff file or buffer size.

To create the bed files that will contain the non-genic regions to be excluded for each reference, use the script ./scripts/code/make_nongenic_bed.py:

    python make_nongenic_bed.py -gff <gff_file> -b <buffer_size> > <reference_name>.NonGenicMask.bed

The -b option specifies a buffer on either side of gene boundary.  E.g Gene1 ends at 100bp and Gene2 starts at 500bp.  If b=10, then the non-genic region will be 110-490.  Current analysis are being done with -b 2000.

While Lumpy requires a bed file specifying regions to *exclude*, Genome STRiP wants a _.list_ file of regions to *include*.  To make these _.list_ files, do:

    python make_genic_bed.py -gff <gff_file> -b <buffer_size> > <reference_name>.genic.list

### Lumpy

The team behind Lumpy has created several software packages (speedseq, svtyper, and svtools) all of which are an encouraged or necessary part of a pipeline to go from read data to a final multi-sample VCF.

We implement Lumpy on a per-individual basis via the speedseq sv command, which will subsequently run svtyper and CNVnator, including all results in a single output file per individual.  Subsequent steps (implemented ) latter steps are important for merging variant calls across individuals using svtools.

NOTE: See instructions above (Requirements) for installing ROOT prior to proceeding

Prior to running speedseq sv, you must also split each fasta by chromosome and label them exactly as the sequence name appears in the original fasta (filenames must end with ".fa").  For each reference, create a directory containing ONLY these chromosome fastas.  Make one copy of the ./speedseq/bin/speedseq.config file per reference genome and provide the directory containing the chromosome fastas after CNVNATOR_CHROMS_DIR=.  Name these config files speedseq.REF.config where the REF name matches that used in the reference_path_file.

Note: Make certain to run speedseq sv using the following options: -v -d -P -g -k option as subsequent steps will utilize CNVnator files in the temporary directories, assume that SVTyper has been run and require LUMPY's probability curves.

To generate commands for speedseq sv, use the script ./scripts/code/Generate_SpeedSeqSV_commands.py:

    python Generate_SpeedSeqSV_commands.py -b <bam_directory> \
      -r <reference_path_file> \
      -o <output_directory> \
      -c <number_of_cores> \
      -s <speedseq_directory> > <command_file>

The reference_path_file is the same as before, but with an additional column containing the path to the non-genic bed files created with make_nongenic_bed.py.

Once all files have run through speedseq sv, we can start the process of merging variants across individuals to create the final multi-sample VCF.  These steps will be performed using a number of functions within svtools.  svtools is single-threaded so parallelization is implemented via GNU parallel in all corresponding shell scripts.

You can download and install svtools plus all requirements with:

    source activate py27
    conda install -c bioconda svtools

NOTE: You must activate this virtual environment every time you wish to run svtools.

After everything installs successfully, the first step is to use svtools lsort, which will combine all of the individual VCF files produced from speedseq sv into a single, sorted VCF followed by svtools lmerge, which will merge overlapping SV calls.

These two processes are implemented in ./scripts/jobs/SVtools_SortAndMerge.sh as follows:

    qsub SVtools_SortAndMerge.sh -F "<input_vcf_dir> <output_vcf> <temp_dir> <batch_size> <percent_slop>” \
    -l mem=20gb,walltime=24:00:00 \
    -q mesabi \
    -o <stdout_file> \
    -e <stderr_file>

All arguments are required, must follow correct order, and must be space-delimited within quotations following the -F flag.  You can change the requested memory and walltime using the -l flag.  

_input_vcf_dir_ is directory containing per-individual VCFs.  If you have called SVs against multiple reference genomes, create separate directories for each reference and call the script once per input directory.  output_vcf is the full path to output filed cannot already exist.  temp_dir is where to place temporary files during the sorting step.  batch_size is the number of vcfs to sort at a time, which determines amount of memory required.  percent_slop is the percent increase of the breakpoint confidence interval both up and down stream, which determines the aggressiveness of the merging process. 

The next step is to run the svtools genotype and svtools copynumber.  Prior to this, you must create the coordinates_file needed for CNVnator for each of the merged vcf files.

    zcat <merged_ref.vcf.gz> > <merged_ref.vcf>
    source activate py27
    create_coordinates -i <merged_ref.vcf> -o <ref_coordinates_file>
_create_coordinates_ is a script that ships with svtools and should install as an executable in your /usr/bin/ directory.  It will be available after activating the python 2.7 virtual environment in which svtools was installed.

To generate commands for svtools genotype and svtools copynumber use the script in ./scripts/code/Generate_SVtools_Genotype_commands.py

    python Generate_SVtools_Genotype_commands.py \
      -b <bam_directory_used_in_speedseq_sv> \
      -o <output_directory> \
      -v <merged_ref.vcf.gz> \
      -c <ref_coordinates_file> \
      -s <path_to_speedseq_software_directory> \
      -w <window_size> >> <svtools_genotype_commands_file>
The default value of 300 for the window size parameter seems to perform robustly.  Lower values will likely result in job failures.  You will need to run this once for each merged_ref.vcf.  These commands can be run using ./scripts/jobs/SVtools_Genotype.sh, which uses GNU parallel.  They can also be split and run as multiple GNU parallel jobs via a task array as we did with Speedseq_large.sh.

NOTE: The ‘>>’ will append commands to a file, so if an older file already exists, you should delete the older file before running this script.

The final step in the Lumpy/svtools pipeline involves pasting the individual genotyped vcfs that include information from CNVnator, further prune the combined VCF for redundant SV calls, and classify remaining calls as valid or artificial.  Prune takes a clustering approach based on a specified evaluation parameter, and prunes sites in the same cluster that are within a given distance.

Use the script ./scripts/jobs/SVtools_PastePruneClassify.sh:

    qsub SVtools_PastePruneClassify.sh -F "<vcf_dir> <merged.vcf> <out_prefix> <dist> <eval_param> <te.bed>”
    -l mem=20gb,walltime=24:00:00
    -o stdout_file \
    -e stderr_file

vcf_dir is the input directory containing per-individual VCFs.  If you have called SVs against multiple reference genomes, create separate directories for each reference and call the script once per input directory.  merged.vcf is the merged vcf that was used in svtools genotype/copynumber.  out_prefix is the prefix to use for output vcf.  These files cannot already exist.  dist is the max separation distance (bp) of adjacent loci in a cluster [50].  eval_param is the evaluating parameter for choosing best bedpe in a cluster(e.g. af=AlleleFrequency default:af).  Currently, it seems that “af” is the only option for eval_param.  te.bed is a gzipped bed file from repeat masker or elsewhere that specifies te locations.

### Genome STRiP

**NOTE**: Genome STRiP was abandoned, so the following section can be ignored.

Many of the Genome STRiP utilities rely on a script called Queue.jar that will automatically launch and manage a large number of (grand)child processes.  In order for these to run, they must inherit the environment from the parent job.  One way to accomplish this is to modify your _~/.bashrc_ profile, so that all necessary paths and softwares are loaded every time a job is launched under your username.  For me, the following lines were added to _~/.bashrc_:

    module load htslib/1.6
    module load samtools
    module load liblzma
    module load java/jdk1.8.0_144
    module load libdrmaa/1.0.13
    SV_DIR="/home/hirschc1/pmonnaha/software/svtoolkit"
    export LD_LIBRARY_PATH=${SV_DIR}:${LD_LIBRARY_PATH}
    export SV_DIR
    export PATH=${SV_DIR}:${PATH}
    export LD_LIBRARY_PATH=/panfs/roc/msisoft/libdrmaa/1.0.13/lib/:${LD_LIBRARY_PATH}

Prior to running Genome STRiP, you also have to create a “MetaData Bundle” for each reference.  See [Genome STRiP: Preparing a reference](scripts/GSTRiP_RefPrep.md)

Once the MetaData bundle is complete, the CNVDiscoveryPipeline can be run with *GSTRiP_CNVDiscoveryPipeline.sh*. All arguments are required and must 
follow correct order.  User should also specify queue and requested resources via command line.  See note below for use of default values.

Usage:

    qsub GSTRiP_CNVDiscovery.sh -F \"refName outDir memGb bamListFile tilingWindowSize tilingWindowOverlap maximumReferenceGapLength boundaryPrecision minimumRefinedLength\"
refName: e.g. W22, B73, PH207, or PHB47.  Paths are hardcoded.  Change if necessary.
outDir: Output Directory
memGb: number of Gb of RAM for CNVDiscovery to request.  This should be about 2Gb less than requested from the HPC (e.g. qsub ... -l mem=...)
bamListFile: file containing list of bams to include in the analysis.  Must end with _.list_

Remaining parameters will control sensitivity and should be calibrated for the mean coverage across samples.  Generally, larger values will reduce sensitivity and increase the minimum detectable size of SVs (the final parameter listed below).  These correspond to the _tilingWindowSize_, _tilingWindowOverlap_, _maximumReferenceGapLength_, _boundaryPrecision_, and _minimumRefinedLength_, respectively.  See the [CNVDiscoveryPipeline Documentation](http://software.broadinstitute.org/software/genomestrip/org_broadinstitute_sv_qscript_CNVDiscoveryPipeline.html) for more information on the meaning of these parameters.

If you have samples of varying coverage, consider trying a range of different values while subsetting your data to include only samples with sufficient coverage for the current values (see this [post](https://gatkforums.broadinstitute.org/gatk/discussion/12267/window-size-for-samples-of-varying-coverage#latest) for a discussion on this issue). If you simply want to run default settings (for samples at approx. 25x coverage), add these values at the end of the -F command:

    1000 500 1000 100 500

For example:

    qsub ~/JobScripts/GSTRiP_CNVDiscoveryPipeline.sh -l mem=32gb,walltime=96:00:00 -F "W22 /panfs/roc/scratch/pmonnaha/Maize/gstrip/w22_defaults 30 /home/hirschc1/pmonnaha/misc-files/gstrip/W22_E2_Bams.list 1000 500 1000 100 500" -A hirschc1 -q mesabi

For my initial attempts, I ran default settings with all samples included as well as a low-sensitivity run (again with all samples) and a high-sensitivity run, including only samples with approx. >40x coverage. Low-sensitivity parameters (calibrated for 10x coverage) are:

    2000 1000 1500 150 1000

High-sensitivity parameters for >40x coverage are:

    500 250 500 75 250
High-coverage samples used in these runs can be found in _./accessory/High_cov_samples.txt_

#### Filtering redundant calls

Genome STRiP includes a utility (RedundancyAnnotator) that can filter the output of CNVDiscoveryPipeline for redundant calls based on overlapping start and end coordinates.  See the [Redundancy Annotator Documentation](http://software.broadinstitute.org/software/genomestrip/org_broadinstitute_sv_annotation_RedundancyAnnotator.html) for more detail.

For a pair of Genome STRiP VCFs, run the RedundancyAnnotator with the job script _./scripts/jobs/GSTRiP_RedundancyAnnotator.sh_

Usage: 

    qsub GSTRiP_RedundancyAnnotator.sh -F \"VCFile refName outPath memGb duplicateOverlapThreshold VCF2File\"
_VCFile_: vcf that you would like ot genotype the samples in bamListFile
_refName_: e.g. W22, B73, PH207, or PHB47.  Paths are hardcoded.  Change if necessary.
_outPath_: Full path to output file
_memGb_: number of Gb of RAM to request.  This should be about 2Gb less than requested from the HPC (e.g. qsub ... -l mem=...)
_VCF2File_: VCF to compare with original

To iteratively run this shell script on a collection of VCFs (e.g. resulting from different sensitivity settings), use the wrapper script _./scripts/code/GSTRiP_RedundancyAnnotator_wrapper.py_.  This will compare each VCF with itself and all others (parallelized by chromosome), extract the preferred variant according to GSTRiP, and concatenate the results into a single VCF.

Usage: 

    GSTRiP_RedundancyAnnotator_wrapper.py -f vcf_file_list \
        -p path_to_sh \
        -m memGB \
        -r refName \
        -d duplicateOverlapThreshold \
        -o output_VCF_name \
        -v vcftools_perl_folder \
-f : full path to each VCF to compare.  One path per line.
-p : Path to GSTRiP_RedundancyAnnotator.sh
-m : Amount of memory to request from HPC
-r : W22, B73, PH207, or PHB47.  Paths are hardcoded.  Change if necessary.
-d : Minimum overlap required to determine two events are redundant
-o : output file name
-v : path to directory containing the vcftools perl modules

*NOTE*: If you are running this wrapper for PH207 as reference, you need to add the flag _--chromosome10_ in order to accommodate chromosome naming convention for this reference.

## SV Consolidation

The purpose of the code described in this section is to match SVs in VCFs called against different reference genomes, based on associating homologous genes that have been annotated for each variant.  At this stage, the user should have multiple VCFs corresponding to the SAME samples mapped to different reference genomes.  If the same samples are not present across VCFs, the following steps will fail.  Additionally, it is assumed that the user will have run the Gene-key pipeline (https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline; specifically, the output of "All_By_All_compare"), which links homologous genes across maize annotations.  

The main task of linking SVs across reference genomes is done by the script _scripts/code/SVMap.py_, but there are at least two necessary steps to prepare (i.e. format and compress for fast access) the two major components mentioned above (the collection of VCFs, and the collection of gene-key files).  Use the _annotate_pickle_vcf.py_ script in _scripts/code_ to prepare the collection of VCFs.  

    usage: annotate_pickle_vcf.py [-h] -f vcf_info_file -o output_directory -s
                              output_suffix [-sp SURVIVOR_ant_path]
                                [-ad annotation_distance] [-b buffer] \\
    This script prepares the data contained in multiple SV VCFs (corresponding to
    the same samples called against different reference genomes) for subsequent
    analysis with SVMap.py (which links the SVs across reference genomes).
    optional arguments:
      -h, --help            show this help message and exit
      -f vcf_info_file      tab delimited file with each line containing 1.) the
                            reference genotype ID (e.g. B73; used for naming), 2.)
                            a bed file with gene locations ONLY and 3.) the vcf
                            file.
      -o output_directory   Output Directory
      -s output_suffix      Output Suffix
      -sp SURVIVOR_ant_path
      -ad annotation_distance
                            distance from SV to buffer for looking for gene
                            overlap when annotating merged vcf with gene info;
                            this is accomodated by -b, which provides more
                            explicit info regarding the location where the overlap
                            is found, so this can be left at 0 (default)
See _vcf_list_AnnotatePickle.txt_ in _accessories_ for an example of the `<vcf_info_file>`.

To prepare the collection of gene-key files, use the script code/pickle_homologue_dicts.py

    usage: pickle_homologue_dicts.py [-h] -f geneKey_fileList -o output_directory \\
    This script prepares the gene-keys produced by JM's gene-key pipeline (https:/
    /github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline), so
    that they can be used to link SVs across references (SVMap.py) \\ \\ 
    optional arguments:
      -h, --help           show this help message and exit
      -f geneKey_fileList  tab delimited file 3 columns: reference1 reference2
                       ref1-ref2_gene-key-file. Note that there is a
                       directionality to the gene-key files. E.g. The gene-key
                       file, 500kb_W22_B73_AllbyAll_res.txt, specifically
                       corresponds to W22 as reference1 and B73 as reference2
      -o output_directory
See scripts/accessory/geneKey_fileList.txt for the list of gene-key files that were used.

### SVMap

This is the primary script that is used to cross-validate SV variants across
reference genomes. Users should have run the necessary scripts to prepare the
VCFs (output of annotate_pickle_vcfs.py) and gene-keys (output of
pickle_homologue_dicts.py). The first output file (suffix: .results.txt)
summarizes relavent info for each match, whereas the second output file
(suffix: .samples.txt) contains the list of samples that were filtered and/or
did not match across references for each match, which could be used for
filtering sample genotypes when creating consensus calls across references. 

    usage: SVMap.py [-h] -f npz_file -o output_directory [-s output_prefix]
                [-gq genotype_quality_threshold] [-mi minimum_individuals]
                [-d tmp_directory] \\
    optional arguments:
      -h, --help            show this help message and exit
      -f npz_file           tab delimited file with 3 colums: 1.)key name used in
                            homologue dictionary (e.g. B73, W22, PHB47, and PH207)fe
                            2.) the npz file of the converted vcf (created with
                            annotate_pickle_vcf.py), and 3.) the compressed
                            dictionary containing homolog info (created by
                            pickle_homologue_dicts.py)
      -o output_directory   Directory to write the two output files
      -s output_prefix
      -gq genotype_quality_threshold
                            Filter all genotypes below this threshold. Directly
                            impacts matching and subsequently genotype distance
                            calculations; Default = 0
      -mi minimum_individuals
                            minimum number of inidividuals (with sufficient
                            genotype quality in both references) needed to
                            calculate a distance between genotype matrices;
                            Default = 0
      -d tmp_directory

## Notes

* * *

### Parallelization

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

### Useful Commands

Get deduplicate percentages from fastqc

    grep Deduplicate <fastqc_results_dir>/*/fastqc_data.txt | awk '{print $4}' > <output_file>


Get adaptor info from cut adapt

    grep -n -e filter -e Read <stdout_directory>/cutadapt_set*.o-* | grep "%" | awk '{print $1, $2, $NF}' | cut -d ":" -f 1,3 | sed 's/written//' | sed  's/(//' | sed 's/%)//' | sed 's/\:/\t/' | rev | cut -d "_" -f 1 | rev > <output_file>


Get sickle results

    grep -n -e kept -e discarded <stdout_directory>/cutadapt_set*.o* | cut -d ":" -f 1,3,4 | awk '{print $1, $2, $4, $5}' | rev | cut -d "_" -f 1 | rev > <output_file>


Get percent of duplicates based on mapping.  The error files from the SpeedSeq align commands contain a line that has the percent of duplicates for the mapped fastq.  To retrieve these values along with sample names use:

    grep duplicates Speedseq_large.* | awk '{split($0,a,"CutAdapt/") ; print(a[2],a[3],a[4])}' | cut -d ")" -f 1 -s | awk '{split($1,a,"_R1");print(a[1],$NF)}' | cut -d "%" -f 1 -s | sed 's/(/ /'


#### Merging newly-mapped BAMs with pre-existing, already-merged BAMs.  

Run the below python script to generate the commands to be subsequently run on the cluster.  

    python scripts/code/Generate_MergeExistingBAMs_commands.py \
        -u <unmerged_bam_directory> \
        -m <merged_bam_directory> \
        -k <merge_Key> \
        -o <output_bam_directory> \
        -r <Reference_Path_Key> \
        -s <path_to_speedseq_directory> > <command_file>

The *merge_key* is essentially the same as the *sample_fastq_key* used previously, except I subset it to include only the necessary sample IDs.  You can use the shell script *./scripts/jobs/Sambamba_MergeBAMs.sh* to submit the jobs in the *command_file* as a task array, via:

    qsub -t <XX>-<YY> Sambamba_MergeBAMs.sh -F "<path_to_command_file>"

, where XX and YY correspond to the range of line numbers in the *command_file* that user has specified.


  

