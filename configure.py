#! /usr/bin/env python
#-----------------------------------------------------------------------------------------
import argparse
import glob
import re

# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

# Step 1. Prepare parser, add each of the arguments
parser = argparse.ArgumentParser()

user_directory = 'user/'
user_choices = glob.glob(user_directory + '*.F90')
user_choices = [choice[len(user_directory):-4] for choice in user_choices]
parser.add_argument('--user',
    default='user_weibel',
    choices=user_choices,
    help='select user file')

parser.add_argument('-hdf5',
    action='store_true',
    default=False,
    help='enable HDF5 & use h5pfc compiler')

parser.add_argument('-mpi',
    action='store_true',
    default=False,
    help='enable mpi & use mpif90 compiler')

parser.add_argument('-intel',
    action='store_true',
    default=False,
    help='enable intel compiler')

parser.add_argument('-debug',
    action='store_true',
    default=False,
    help='enable DEBUG flag')

parser.add_argument('-3d',
    action='store_true',
    default=False,
    help='enable 3d')

args = vars(parser.parse_args())

# Step 2. Set definitions and Makefile options based on above arguments

makefile_options = {}
makefile_options['USER_FILE'] = args['user']

makefile_options['COMPILER_COMMAND'] = ''
makefile_options['COMPILER_FLAGS'] = ''
makefile_options['PREPROCESSOR_FLAGS'] = ''

if args['hdf5']:
    makefile_options['COMPILER_COMMAND'] += 'h5pfc '
    makefile_options['PREPROCESSOR_FLAGS'] += '-DHDF5 '
elif args['mpi']:
    makefile_options['COMPILER_COMMAND'] += 'mpif90 '
    makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI '
else:
    makefile_options['COMPILER_COMMAND'] += 'gfortran '

if args['debug']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG -fcheck=all -fimplicit-none -fbacktrace '

if args['debug'] and args['intel']:
    makefile_options['COMPILER_FLAGS'] += '-O3 -qopenmp-simd -qopt-report=5 '

if args['intel']:
    makefile_options['MODULE'] = '-module '
else:
    makefile_options['MODULE'] = '-J '

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
