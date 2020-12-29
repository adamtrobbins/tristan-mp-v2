#!/usr/bin/env python3
import sys
import os
import glob
import shutil
from abc import ABC, abstractmethod

from datetime import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', action='store_true', default=False, help='only test compilation.')
parser.add_argument('-v', action='store_true', default=False, help='verbose mode (full output).')
parser.add_argument('-d', action='store_true', default=False, help='diagnostic mode.')
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
class Simulation(ABC):
  path = None
  exe = None
  exe_full = None
  submit = None
  submit_full = None
  input = None
  input_full = None
  def __init__(self, flags, walltime='00:30:00', nproc=28, params={}):
    self.walltime = walltime
    self.flags = flags
    self.params = params
    self.nproc = nproc
  @abstractmethod
  def diag(self, ax):
    pass

class TwoStream(Simulation):
  jobid = 'twostream'
  userfile = 'user_twostream'
  dimension = 1
  def diag(self, ax):
    hist = isolde.parseHistory(self.path + '/output/history')
    omegap0 = self.params['algorithm']['c'] / self.params['plasma']['c_omp']
    rate = 0.5 * (0.5)**0.5 / (self.params['problem']['shift_gamma'])**(1.5)
    time = hist['time'] * omegap0; E2 = hist['E^2'] / hist['Etot'][0]; ax.plot(time, E2)
    xs = np.linspace(10, 40, 10); ys = np.exp(2 * rate * xs); ys = (E2[time>10][0]) * (ys / ys[0]); ax.plot(xs, ys)
    ax.set_ylim(1e-4, 1e-1); ax.set_xlim(0, 200); ax.set_yscale('log');
    ax.set_xlabel(r'$t\omega_{\rm p0}$'); ax.set_ylabel(r'$U_E / E_{\rm tot}$')
    ax.axvline(ax.get_xlim()[0], color='black'); ax.axhline(ax.get_ylim()[0], color='black')
    ax.set_title(self.jobid)

class PlasmaOsc(Simulation):
  jobid = 'plasmaosc'
  userfile = 'user_langmuir'
  dimension = 1
  def diag(self, ax):
    exs = []; steps = np.arange(50)
    for step in steps:
      flds = isolde.getFields(self.path + '/output/flds.tot.%05d' % step)
      exs.append(flds['ex'][0, 0, int(self.params['grid']['mx0'] / 4)] / 1e-4)
    ax.plot(steps * self.params['output']['interval'] / (2 * np.pi * self.params['plasma']['c_omp'] / 0.45), exs)
    ax.set_xlabel(r'$t\omega_{\rm p0}$'); ax.set_ylabel(r'$E_x$')
    ax.set_title(self.jobid)

# Here specify the test simulations and give additional specs of the environment
common_flags = ' -perseus -hdf5 -debug'
simulations = [
               TwoStream(common_flags,
                          params={
                            'node_configuration': {'sizex' : 28},
                            'time': {'last' : 50000},
                            'grid': {'mx0' : 5600, 'tileX' : 100},
                            'algorithm': {'c' : 0.35, 'nfilter': 8},
                            'output': {'enable' : 0, 'hst_enable' : 1, 'hst_interval' : 5},
                            'plasma': {'ppc0' : 64, 'sigma' : 10, 'c_omp' : 40},
                            'particles': {'nspec' : 2, 'maxptl1' : 1e5, 'm1' : 1, 'ch1' : -1, 'maxptl2' : 1e5, 'm2' : 1, 'ch2' : -1},
                            'problem': {'shift_gamma' : 2.5}}
                          ),
               PlasmaOsc(common_flags, nproc=4,
                          params={
                            'node_configuration': {'sizex' : 4},
                            'time': {'last' : 1000},
                            'grid': {'mx0' : 1120, 'tileX' : 10},
                            'algorithm': {'nfilter': 8},
                            'output': {'enable' : 1, 'interval': 10, 'hst_enable' : 1, 'hst_interval' : 5},
                            'plasma': {'ppc0' : 500, 'sigma' : 1, 'c_omp' : 10},
                            'particles': {'nspec' : 2, 'maxptl1' : 1e8, 'm1' : 1, 'ch1' : -1, 'maxptl2' : 1e8, 'm2' : 1, 'ch2' : 1},
                            'problem': {'upstream_T' : 1e-5, 'amplitude' : 0.01, 'nwaves' : 1}}
                         )]

if (options.d):
  # diagnostic mode where you analize the test results
  import matplotlib.pyplot as plt
  import numpy as np
  import tristanVis.isolde as isolde
  from matplotlib import rc
  plt.style.use('fivethirtyeight')
  rc('font',**{'family':'monospace', 'size':15})
  rc('text', usetex=True)
  fig = plt.figure(figsize=(16, 10))
  for ii, simulation in enumerate(simulations):
    ax = plt.subplot(2, 2, ii + 1)
    simulation.path = testdir_full + '/%02d_' % (ii + 1) + simulation.jobid
    simulation.diag(ax)
  plt.show()
else:
  # regular mode where you compile and run tests
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

          sub.write('DIR={}\n'.format(simulation.path))
          sub.write('EXECUTABLE=$DIR/{}\n'.format(simulation.exe))
          sub.write('INPUT=$DIR/{}\n'.format(simulation.input))
          sub.write('OUTPUT_DIR=$DIR/output\n'.format(simulation.path))
          sub.write('SLICE_DIR=$DIR/slices\n'.format(simulation.path))
          sub.write('REPORT_FILE=$DIR/report\n'.format(simulation.path))
          sub.write('ERROR_FILE=$DIR/error\n\n'.format(simulation.path))

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
