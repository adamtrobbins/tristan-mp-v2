#!/usr/bin/env python3
import sys
import os
import glob


# Snipped specifying the modules used for compilation
modules = ['intel-mkl/2017.4/5/64', 
            'intel/17.0/64/17.0.5.239', 
            'intel-mpi/intel/2017.5/64',
            'hdf5/intel-17.0/intel-mpi/1.10.0']

# Class specifying the environment for each test simulation (tristan & slurm)
class TestSim:
    walltime = '24:00:00'
    ppn = 1
    nodes = 1
    dims = 1
    conf = 'python configure.py -perseus -hdf5 --nghosts=5'
    def __init__(self, name):
            self.name = name


# Here specify the test simulations and give additional specs of the environment
testnames = [TestSim('twostream')]
testnames[0].ppn = 28

workdir = os.getcwd()
outdir = '/scratch/gpfs/mahlmann/Tristan/tristan-sims'

if not os.path.exists(f'{outdir}/test'):
    os.makedirs(f'{outdir}/test')

# Clean up the testing environment from previous tests
os.system(f'rm -rf {outdir}/test/*')

with open(f'{outdir}/test/test.log', 'w+') as tlog:

    for testname in testnames:

        tlog.write(f'{testname.name}'.ljust(30, '.'))

        # Make clean
        os.system('make clean')

        # Configure tristan for the specific simulation
        testname.conf = testname.conf+f' -{testname.dims}d --user user_{testname.name}'
        os.system(testname.conf)

        # Try to compile tristan with the desired configuration
        os.system('make')

        # Check if compilation has been successfull. If yes, write a submitscript
        if os.path.isfile(f'{workdir}/exec/tristan-mp{testname.dims}d'):

            # Write submitscript for that specific simulation
            with open(f'{workdir}/submit_{testname.name}', 'w+') as f:
                f.write(f'#!/bin/bash\n')
                f.write(f'#SBATCH --time={testname.walltime}\n')
                f.write(f'#SBATCH -N {testname.nodes} -n {testname.ppn*testname.nodes}\n')
                f.write(f'#SBATCH --ntasks-per-node={testname.ppn}\n')
                f.write(f'#SBATCH --cpus-per-task=1\n')
                f.write(f'#SBATCH --output={outdir}/test/{testname.name}/tristan-v2.out\n')
                f.write(f'#SBATCH --error={outdir}/test/{testname.name}/tristan-v2.err\n')
                f.write(f'\n')

                for mod in modules:
                    f.write(f'module load {mod}\n')

                f.write(f'srun -N {testname.nodes} -n {testname.ppn*testname.nodes} --ntasks-per-node={testname.ppn} exec/tristan-mp{testname.dims}d -i {workdir}/inputs/input.{testname.name} -o {outdir}/test/{testname.name}')

            os.mkdir(f'{outdir}/test/{testname.name}')
            os.system(f'sbatch {workdir}/submit_{testname.name}')
            os.remove(f'{workdir}/submit_{testname.name}')
            tlog.write('[OK]\n')

