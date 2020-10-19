#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -N 1 -n 4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/gpfs/mahlmann/Tristan/tristan-sims/WaterbagTest_TESTF/tristan-v2.out
#SBATCH --error=/scratch/gpfs/mahlmann/Tristan/tristan-sims/WaterbagTest_TESTF/tristan-v2.err

module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16

srun -N 1 -n 4 --ntasks-per-node=4  exec/tristan-mp2d -i inputs/input.2d_Waterbag -o /scratch/gpfs/mahlmann/Tristan/tristan-sims/WaterbagTest_TESTF
