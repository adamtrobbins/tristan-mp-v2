#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH -N 4 -n 112
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/gpfs/mahlmann/Tristan/tristan-sims/REC_FASTTWO_2_CRAND/tristan-v2.out
#SBATCH --error=/scratch/gpfs/mahlmann/Tristan/tristan-sims/REC_FASTTWO_2_CRAND/tristan-v2.err

module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16

time srun -N 4 -n 112 --ntasks-per-node=28  exec/tristan-mp2d -i inputs/input.2d_rec -o /scratch/gpfs/mahlmann/Tristan/tristan-sims/REC_FASTTWO_2_CRAND
