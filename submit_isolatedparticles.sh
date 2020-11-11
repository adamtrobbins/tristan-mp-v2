#!/bin/sh
#SBATCH --time=01:00:00
#SBATCH -N 1 -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/gpfs/mahlmann/Tristan/tristan-sims/Isoltaed_Test_CCOM/tristan-v2.out
#SBATCH --error=/scratch/gpfs/mahlmann/Tristan/tristan-sims/Isoltaed_Test_CCOM/tristan-v2.err

module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16

time srun -N 1 -n 1 --ntasks-per-node=1  exec/tristan-mp2d -i inputs/input.2disolatedparticles -o /scratch/gpfs/mahlmann/Tristan/tristan-sims/Isoltaed_Test_CCOM
