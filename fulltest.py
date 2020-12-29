#!/usr/bin/env python3
import sys
import os
import glob
import shutil
from datetime import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', action='store_true', default=False, help='only test compilation.')
parser.add_argument('-v', action='store_true', default=False, help='verbose mode (full output).')
options = parser.parse_args()

suffix = '' if options.v else ' >/dev/null 2>&1'

# modules for compilation
modules = ['intel-mkl/2017.4/5/64',
            'intel/17.0/64/17.0.5.239',
            'intel-mpi/intel/2017.5/64',
            'hdf5/intel-17.0/intel-mpi/1.10.0']

# global variables
outdir = "/scratch/gpfs/hakobyan/tristan_2/test"

codedir = os.getcwd()
testdir = 'test_0'
# testdir = 'test_' + datetime.now().strftime("%H.%M_%d.%m.%Y")
testdir_full = outdir + '/' + testdir

# simulation environment (tristan & slurm)
class Simulation:
  path = None
  exe = None
  exe_full = None
  submit = None
  submit_full = None
  input = None
  input_full = None
  def __init__(self, jobid, userfile, dimension, flags,
               walltime='00:30:00', nproc=28, params={}):
    self.jobid = jobid
    self.walltime = walltime
    self.nproc = nproc
    self.dimension = dimension
    self.userfile = userfile
    self.flags = flags
    self.params = params

# Here specify the test simulations and give additional specs of the environment
common_flags = ' -perseus -hdf5 -debug'
simulations = [Simulation('twostream', 'user_twostream', 1, common_flags,
                          params={
                              'node_configuration': {'sizex' : 28},
                              'time': {'last' : 50000},
                              'grid': {'mx0' : 5600, 'tileX' : 100},
                              'algorithm': {'c' : 0.35, 'nfilter': 8},
                              'output': {'enable' : 0, 'hst_enable' : 1, 'hst_interval' : 5},
                              'plasma': {'ppc0' : 64, 'sigma' : 10, 'c_omp' : 40},
                              'particles': {'nspec' : 2, 'maxptl1' : 1e5, 'm1' : 1, 'ch1' : -1, 'maxptl2' : 1e5, 'm2' : 1, 'ch2' : -1},
                              'problem': {'shift_gamma' : 2.5}
                            })]

if not os.path.exists(testdir_full):
  os.makedirs(testdir_full)

with open(testdir_full + '/test.log', 'w+') as testlog:
  # load modules
  for ii, simulation in enumerate(simulations):
    # create directory for simulation
    simulation.path = testdir_full + '/%02d_' % (ii + 1) + simulation.jobid
    if os.path.exists(simulation.path):
      shutil.rmtree(simulation.path)
    os.makedirs(simulation.path)

    testlog.write(('TEST_#{}_'.format(ii+1) + simulation.jobid).ljust(50, '.') + '\n')

    # configure
    config_command = 'python configure.py '
    config_command += simulation.flags
    config_command += ' -{}d'.format(simulation.dimension)
    config_command += ' --user=' + simulation.userfile
    os.system(config_command + suffix)

    # clean
    os.system('make clean' + suffix)
    # compile
    os.system('make all' + suffix)

    # check if compilation successfull
    simulation.exe = 'tristan-mp{}d'.format(simulation.dimension)
    simulation.exe_full = simulation.path + '/' + simulation.exe
    if os.path.isfile(codedir + '/exec/' + simulation.exe):
      testlog.write('compilation'.ljust(46, '.') + '[OK]\n')
      print ('Compilation of `{}` done.'.format(simulation.jobid))

      # move executable
      os.system('mv {} {}'.format(codedir + '/exec/' + simulation.exe, simulation.path) + suffix)

      # clean
      os.system('make clean' + suffix)

      # write input
      simulation.input = 'input.' + simulation.jobid
      simulation.input_full = simulation.path + '/' + simulation.input
      with open(simulation.input_full, 'w+') as inp:
        for block in simulation.params.keys():
          inp.write('\n<{}>\n\n'.format(block))
          for var in simulation.params[block].keys():
            inp.write('  {}  =  {}\n'.format(var, simulation.params[block][var]))
      testlog.write('input file'.ljust(46, '.') + '[OK]\n')

      # write submit
      simulation.submit = 'submit_' + simulation.jobid
      simulation.submit_full = simulation.path + '/' + simulation.submit
      with open(simulation.submit_full, 'w+') as sub:
        sub.write('#!/bin/bash\n')
        sub.write('#SBATCH -t {}\n'.format(simulation.walltime))
        sub.write('#SBATCH -n {}\n'.format(simulation.nproc))
        sub.write('#SBATCH -J {}\n'.format(simulation.jobid))
        sub.write('#SBATCH --output={}/tristan-v2.out\n'.format(simulation.path))
        sub.write('#SBATCH --error={}/tristan-v2.err\n\n'.format(simulation.path))

        sub.write('EXECUTABLE={}\n'.format(simulation.exe_full))
        sub.write('INPUT={}\n'.format(simulation.input_full))
        sub.write('OUTPUT_DIR={}/output\n'.format(simulation.path))
        sub.write('SLICE_DIR={}/slices\n'.format(simulation.path))
        sub.write('REPORT_FILE={}/report\n'.format(simulation.path))
        sub.write('ERROR_FILE={}/error\n\n'.format(simulation.path))

        for module in modules:
          sub.write('module load {}\n'.format(module))

        sub.write('\nmkdir $OUTPUT_DIR\n\n')

        sub.write('srun $EXECUTABLE -i $INPUT -o $OUTPUT_DIR -s $SLICE_DIR -r $RESTART_DIR -R $RESTART > $REPORT_FILE 2> $ERROR_FILE')
      testlog.write('submit file'.ljust(46, '.') + '[OK]\n')
    testlog.write('\n')

  testlog.write('\n\n')
  if not options.c:
    for simulation in simulations:
      os.system('sbatch ' + simulation.submit_full)
      testlog.write(('`{}`'.format(simulation.jobid)).ljust(41, '.') + 'submitted\n')
