#/bin/bash
# ----------- Parameters --------- #
#$ -S /bin/bash
#$ -l mres=5G,h_data=5G,h_vmem=5G
#$ -q lThC.q@compute-10-*
#$ -pe mthread 4
#$ -cwd
#$ -j y
#$ -N {name}
#
# -------- User Variables -------- #
scrtop=/pool/sao/klee
scrdir=$scrtop/$JOB_ID
outputdir=$PBS_O_WORKDIR

module load intel/2018
module load tools/local
module load g16-A

export GAUSS_SCRDIR=$scrdir

#
# ------------- Job -------------- #
#
mkdir $scrdir

echo `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME

g16 < {name}.com > {name}.log

