#!/bin/sh
#SBATCH --time=02:00:00
#SBATCH -N 1 -n 5
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/gpfs/mahlmann/Tristan/tristan-sims/TWOSTREAM_FASTTWO_2_CCOM/tristan-v2.out
#SBATCH --error=/scratch/gpfs/mahlmann/Tristan/tristan-sims/TWOSTREAM_FASTTWO_2_CCOM/tristan-v2.err

module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16

time srun -N 1 -n 5 --ntasks-per-node=5  exec/tristan-mp2d -i inputs/input.1dtwostream -o /scratch/gpfs/mahlmann/Tristan/tristan-sims/TWOSTREAM_FASTTWO_2_CCOM
