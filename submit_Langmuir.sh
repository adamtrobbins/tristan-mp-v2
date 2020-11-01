#!/bin/sh
#SBATCH --time=05:00:00
#SBATCH -N 1 -n 28
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/gpfs/mahlmann/Tristan/tristan-sims/LANGMUIR_FASTTWO_2_NO/tristan-v2.out
#SBATCH --error=/scratch/gpfs/mahlmann/Tristan/tristan-sims/LANGMUIR_FASTTWO_2_NO/tristan-v2.err

module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16

time srun -N 1 -n 28 --ntasks-per-node=28  exec/tristan-mp2d -i inputs/input.langmuir -o /scratch/gpfs/mahlmann/Tristan/tristan-sims/LANGMUIR_FASTTWO_2_NO
