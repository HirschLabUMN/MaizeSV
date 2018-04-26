#!/bin/bash
#PBS -l mem=3gb,nodes=1:ppn=1,walltime=12:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi
#PBS -o /home/hirschc1/pmonnaha/OandE/cutadapt.o
#PBS -e /home/hirschc1/pmonnaha/OandE/cutadapt.e

# Load modules
module load cutadapt

# And the command list
CMD_LIST="/panfs/roc/groups/14/hirschc1/pmonnaha/JobScripts/accessory/Cutadapt_jobs.txt"
SUCCESS_LIST="/panfs/roc/groups/14/hirschc1/pmonnaha/misc-files/Cutadapt_jobs_done.txt"
# taken from the MSI parallel page at
# https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel
# export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"

# parallel --jobs 2 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${CMD_LIST}

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST} | rev | cut -d " " -f 1 | rev >> "${SUCCESS_LIST}"