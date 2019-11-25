#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=18,walltime=30:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -o /home/hirschc1/pmonnaha/OandE/speedseq.o
#PBS -e /home/hirschc1/pmonnaha/OandE/speedseq.e

# Load modules
module load zlib/1.2.8
module load xz-utils/5.2.2_intel2015update3

# And the command list
SCRIPT=`basename $0`
if [ "$#" -ne 1 ]; then
cat << EOF
Usage: sh ${SCRIPT} Command_File
        This script is used to submit speedseq jobs, as listed in Command_File, as a task array.
        
EOF
exit 1
fi

if ! [ -e "$1" ]; then
    echo "Command file:$1 not found" >&2
    exit 1
fi

CMD_LIST="$1"
SUCCESS_LIST="/panfs/roc/groups/14/hirschc1/pmonnaha/misc-files/Successful_Mappings.txt"
# taken from the MSI parallel page at
# https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel
# export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"

# parallel --jobs 2 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${CMD_LIST}

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST} | rev | cut -d " " -f 1 | rev >> "${SUCCESS_LIST}"
