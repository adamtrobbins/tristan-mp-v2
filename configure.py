#! /usr/bin/env python
#-----------------------------------------------------------------------------------------
# configure.py: Athena++ configuration script in python. Original version by CJW.
#
# When configure.py is run, it uses the command line options and default settings to
# create custom versions of the files Makefile and src/defs.hpp from the template files
# Makefile.in and src/defs.hpp.in respectively.
#
# The following options are implememted:
#   -h  --help        help message
#   --prob=name       use src/pgen/name.cpp as the problem generator
#   --coord=xxx       use xxx as the coordinate system
#   --eos=xxx         use xxx as the equation of state
#   --flux=xxx        use xxx as the Riemann solver
#   --nghost=xxx      set NGHOST=xxx
#   -b                enable magnetic fields
#   -s                enable special relativity
#   -g                enable general relativity
#   -t                enable interface frame transformations for GR
#   -debug            enable debug flags (-g -O0); override other compiler options
#   -float            enable single precision (default is double)
#   -mpi              enable parallelization with MPI
#   -omp              enable parallelization with OpenMP
#   -hdf5             enable HDF5 output (requires the HDF5 library)
#   --hdf5_path=path  path to HDF5 libraries (requires the HDF5 library)
#   -fft              enable FFT (requires the FFTW library)
#   --fftw_path=path  path to FFTW libraries (requires the FFTW library)
#   --grav=xxx        use xxx as the self-gravity solver
#   --cxx=xxx         use xxx as the C++ compiler
#   --ccmd=name       use name as the command to call the C++ compiler
#   --include=path    use -Ipath when compiling
#   --lib=path        use -Lpath when linking
#-----------------------------------------------------------------------------------------

import argparse
import glob
import re

# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

# Step 1. Prepare parser, add each of the arguments
parser = argparse.ArgumentParser()

pgen_directory = 'user/'
pgen_choices = glob.glob(pgen_directory + '*.F90')
pgen_choices = [choice[len(pgen_directory):-4] for choice in pgen_choices]
parser.add_argument('--user',
    default='user_weibel',
    choices=pgen_choices,
    help='select user file')

parser.add_argument('-hdf5',
    action='store_true',
    default=False,
    help='enable HDF5 & user h5pfc compiler')

parser.add_argument('-debug',
    action='store_true',
    default=False,
    help='enable DEBUG flag')

parser.add_argument('-3d',
    action='store_true',
    default=False,
    help='enable 3d')

parser.add_argument('-dprec',
    action='store_true',
    default=False,
    help='enable double precision')

args = vars(parser.parse_args())

# Step 2. Set definitions and Makefile options based on above arguments

makefile_options = {}
makefile_options['USER_FILE'] = args['user']

makefile_options['COMPILER_COMMAND'] = 'gfortran ' if (not args['hdf5']) else 'h5pfc '
makefile_options['COMPILER_FLAGS'] = ''
makefile_options['PREPROCESSOR_FLAGS'] = ''

if args['debug']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG '

if args['debug']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG '

if args['dprec']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDPREC '

if args['3d']:
    makefile_options['EXE_NAME'] = 'tristan-mp3d'
    makefile_options['PREPROCESSOR_FLAGS'] += '-DthreeD '
else:
    makefile_options['EXE_NAME'] = 'tristan-mp2d'

# Step 3. Create new files, finish up
with open(makefile_input, 'r') as current_file:
  makefile_template = current_file.read()
for key,val in makefile_options.items():
  makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)
with open(makefile_output, 'w') as current_file:
  current_file.write(makefile_template)

# Finish with diagnostic output
print('Your TRISTAN distribution has now been configured with the following options:')
print('  Userfile:                ' + args['user'])
print('  Dim:                     ' + ('3D' if args['3d'] else '2D'))
print('  Debug mode:              ' + ('ON' if args['debug'] else 'OFF'))
print('  HDF5 output:             ' + ('ON' if args['hdf5'] else 'OFF'))
print('  Compilation command:     ' + makefile_options['COMPILER_COMMAND'] \
    + makefile_options['PREPROCESSOR_FLAGS'] + makefile_options['COMPILER_FLAGS'])
