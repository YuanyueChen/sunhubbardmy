#!/bin/bash

NP=4
ND=1

WORKDIR="$PWD"
echo $WORKDIR
cd $WORKDIR

cat>srun.sub<<endsub
#!/bin/bash
#SBATCH -J srun
#SBATCH -p debug64c512g
###SBATCH --exclusive
#SBATCH -n $NP
#SBATCH --ntasks-per-node=$NP
#SBATCH -t 1:00:00
#SBATCH -o job.%j.%N.out
#SBATCH -e job.%j.%N.err

module load intel-oneapi-compilers
module load intel-oneapi-mkl
module load intel-oneapi-mpi

source run_PC.sh
endsub

sbatch srun.sub