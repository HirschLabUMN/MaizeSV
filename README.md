# WiDiv Structural Variation

#### Summary

Structural variant discovery using whole genome short-read sequence data on 500 inbred lines from the Wisconsin Diversity Panel (WiDiv).  For robust structural variant discovery, we take the consensus variant calls across:


  1. Three state-of-the-art SV discovery tools
       * Lumpy 

       * MetaSV

       * Genome STRiP

  2. At least 4 high-quality reference genomes

       * B73, W22, PH207, and PHB47  


## Fasta preparation

Several subsequent steps in the pipeline work more smoothly if all extraneous scaffolds are removed from the reference sequences.  Use **filter_fasta.py** to achieve this:

To install all necessary components to run **filter_fasta.py**, perform the following steps (taken from <https://www.biostars.org/p/157811/>)

1.Ensure you have python and biopython installed. Type in your terminal:  

    python -c "Import Bio"
    echo $?

If "0" displayed, pass to step 2. In case error message appears, install it by easy_install or pip: 

    sudo easy_install -f <http://biopython.org/DIST/> biopython

or 

    sudo pip install biopython

In case of not having pip installed, type: 

    sudo apt-get install pip


3.Run the program as follows: 

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

Sequence data comes from 3 sources:  

  1. JGI: An initial subset of ~60 lines (see ./misc/JGI_fastqs.txt)

  2. Mark Mikel (MM): High depth sequence data corresponding to 8 lines (MM_fastqs.txt)

  3. Novagene: All remaining samples (Novagene_fastqs.txt)
<br />
  
The sequence data from MM (8 samples listed below) exists in multiple files corresponding to the multiple lanes on which each sample was sequenced.  These files need to be concatenated prior to mapping.

  * B73, G84, G35, J40, G39, LH82, G47, and PH207
<br />

#### Quality Assessment — _Fastqc_

To generate _fastqc_ commands:

    find <fastqdir> -name "*fastq*" -o -name "*fq*" -print | xargs -I {} echo "fastqc --noextract -t <number_of_threads> -o <outdir>" {} > <command_file>

The terms in <> need to be replaced by the user.  

  
You can use **fastqc.sh** in _/scripts/jobs/_ to run these commands in parallel via a task array.  If you are unfamiliar with the different options for parallelizing job submission, see the discussion below in **Notes**.  You will also need to open **fastqc.sh** and modify paths where indicated.

  


#### Adapter and Quality Trimming — _CutAdapt_ and _Sickle_

We use _cutadapt_ and _sickle_ to trim adapters and low quality bases, respectively, from reads.  Both of these softwares need to be downloaded and installed.  Although _cutadapt_ is available as load-able module on MSI, the correct version is not available.

Install _cutadapt_ using pip

    pip install --user --upgrade cutadapt

Install _sickle_ by cloning the git repository and typing ‘make’ from within the directory

    git clone <https://github.com/najoshi/sickle.git>
    cd sickle
    make

To generate commands for cutadapt and sickle, use the script _./scripts/code/_**Generate_CutAdapt_commands.py**.  The use of this script differs depending on whether the paired-end fastq data is in an interleaved format (both forward and reverse reads are interleaved in a single file) versus in separate files.  For interleaved data, simply do:

    python Generate_CutAdapt_commands.py -c <path_to_cutadapt> -s <path_to_sickle> -r <read_path_file> -o <output_directory> >> <commands_file>

whereas for data in separate files do:

python Generate_CutAdapt_commands.py -c path_to_cutadapt -s path_to_sickle -r read1_path_file -r2 read2_path_file -o output_directory >> commands_file

where the read1_path_file contains the forward read files and read2_path_file contains the paths to the reverse reads.  

  


NOTE: the ‘>>’ in the above commands will append the output of each run to the commands_file, which is useful for combining commands for different sets of interleaved/separate fastq files.  However, be sure the commands_file doesn’t already exist, or you may potentially append on top of a bunch of old commands that you don’t want to include.

  


Run commands as a task array with script ./scripts/jobs/Cutadapt.sh.  Open **Cutadapt.sh** and modify the paths where indicated prior to running.

  


Mapping — Speedseq

* * *

We use speedseq for mapping and position-sorting as well as for isolation of split reads and discordant paired reads.  Under the hood, speedseq more efficiently parallelizes bwa mem and pipes the output directly to sambamba and samblaster.  Sambamba is a sam/bam manipulation software with analogous functionality to samtools, like sorting, merging, etc., yet is much faster.  Samblaster efficiently marks duplicate reads and simultaneously extracts discordant and split reads.

  


Install speedseq (also installs _sambamba_ and _samblaster,_ automatically) with:

    module load cmake
    git clone --recursive <https://github.com/hall-lab/speedseq>
    cd speedseq
    make

  


Generate speedseq commands with ./scripts/code/Generate_SpeedSeq_commands.py

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

  


  


 **Merging BAMs — _sambamba_**

* * *

  
The script  **Generate_MergeBAMs_commands.py**  in ./scripts/code/ can be used to generate commands to merge bams with _sambamba:_

python Generate_MergeBAMs_commands.py -k _sample_fastq_key_  \

-o _output_directory_ \

-r reference_path_file \

-c _number_of_cores_  > _sambamba_command_file_

  


Run _sambamba_  commands as task array implemented in ./scripts/jobs/ **Sambamba_MergeBAMs.sh**.  Be sure to open **Sambamba_MergeBAMs.sh** __ and modify paths where necessary.  ALSO, be sure to modify the ‘ppn’ field to reflect the number of cores specified with **Generate_MergeBAMs_commands.py** and also adjust the ‘mem’ field appropriately. 

  


SV Discovery

* * *

Lumpy

  


  


Genome STRiP

  


  


MetaSV

  


  


Notes

* * *

  


Parallelization

There are two options for parallelizing jobs on the MSI cluster:

  1. GNU parallel 

  2. Task Arrays




  


MSI has good help documentation for each approach.  For task arrays, see <https://www.msi.umn.edu/support/faq/how-do-i-use-job-array>.  For GNU parallel, see <https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel>

  


The primary difference is that GNU parallel, for the most part, can only handle single-threaded jobs, whereas task arrays can be used for multi-threaded jobs.  Task arrays also allow for simpler resubmission of failed jobs.

  


Compare the two scripts fastqc_GNUparallel.sh and fastqc_TaskArray.sh to see how job specification differs between the two approaches.  To submit these jobs, you would simply use:

qsub fastqc_GNUparallel.sh

for the former, and:

qsub -t 0-500 fastqc_TaskArray.sh

for the latter, replacing ‘500’ for the total number of commands in your command file.  Failed jobs can be resubmitted with ‘qsub -t <job_number>'

  

