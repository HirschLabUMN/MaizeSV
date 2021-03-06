#!/bin/bash
#PBS -l mem=2976gb,nodes=48:ppn=24,walltime=24:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q large
#PBS -o path/to/stdout/directory/Speedseq_large.o
#PBS -e path/to/stderr/directory/Speedseq_large.e

# Load modules
module load parallel

# MODIFY THE PATH BELOW AS WELL AS -o AND -e PATHS ABOVE
CMD_LIST="/path/to/file/containing/speedseq_commands_${PBS_ARRAYID}"

# taken from the MSI parallel page at
# https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel
export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"

parallel --jobs 1 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${CMD_LIST}
