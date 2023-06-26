#!/bin/bash -l

#Standard output and error:
#SBATCH -o ./out.%j
#SBATCH -e ./err.%j

#Initial working directory:
#SBATCH -D ./

#Job Name
#SBATCH -J pdos-nnso

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=72
#SBATCH --cpus-per-task=1

#Wall clock limit:
#SBATCH --time=0-00:40:00

# Make sure each process only uses one thread, that is, force VASP to use only one thread per core
ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# For pinning threads correctly:
export OMP_PLACES=cores

# Libraries and modules
module load intel/21.7.1
module load impi/2021.7
module load qe/7.2

srun pw.x -i input/scf.in > output/scf.out
srun pw.x -i input/nscf.in > output/nscf.out
srun projwfc.x -i input/projwfc.in > output/projwfc.out