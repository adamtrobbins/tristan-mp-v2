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

unit_directory = 'unit/'
unit_choices = glob.glob(unit_directory + '*.F90')
unit_choices = [choice[len(unit_directory):-4] for choice in unit_choices]

rad_choices = ['no', 'sync', 'ic', 'sync+ic']

# system
parser.add_argument('-perseus',
                    action='store_true',
                    default=False,
                    help='Configure for `Perseus` cluster.')

parser.add_argument('-intel',
                    action='store_true',
                    default=False,
                    help='enable intel compiler')
parser.add_argument('-hdf5',
                    action='store_true',
                    default=False,
                    help='enable HDF5 & use h5pfc compiler')

parser.add_argument('-ifport',
                    action='store_true',
                    default=False,
                    help='enable IFPORT library (`mkdir` etc)')

mpi_group = parser.add_mutually_exclusive_group()
mpi_group.add_argument('-mpi',
                       action='store_true',
                       default=False,
                       help='enable mpi')
mpi_group.add_argument('-mpi08',
                       action='store_true',
                       default=False,
                       help='enable mpi_f08')

# user file
user_group = parser.add_mutually_exclusive_group(required=True)
user_group.add_argument('--user',
                        default=None,
                        choices=user_choices,
                        help='select user file')
user_group.add_argument('--unit',
                        default=None,
                        choices=unit_choices,
                        help='select unit file')

# algorithms
parser.add_argument('--nghosts',
                    action='store',
                    default=3,
                    help='specify the # of ghost cells')

parser.add_argument('-extfields',
                    action='store_true',
                    default=False,
                    help='apply external fields')

parser.add_argument('-absorb',
                    action='store_true',
                    default=False,
                    help='enable absorbing boundaries')

parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable DEBUG flag')

parser.add_argument('-gca',
                    action='store_true',
                    default=False,
                    help='enable GCA mover')

dim_group = parser.add_mutually_exclusive_group(required=True)
dim_group.add_argument('-1d',
                       action='store_true',
                       default=False,
                       help='enable 1d')
dim_group.add_argument('-2d',
                       action='store_true',
                       default=False,
                       help='enable 2d')
dim_group.add_argument('-3d',
                       action='store_true',
                       default=False,
                       help='enable 3d')

parser.add_argument('-dwn',
                    action='store_true',
                    default=False,
                    help='enable particle downsampling')

parser.add_argument('-alb',
                    action='store_true',
                    default=False,
                    help='enable adaptive load balancing')

parser.add_argument('-slb',
                    action='store_true',
                    default=False,
                    help='enable static load balancing')

# extra physics
parser.add_argument('--radiation',
                    default='OFF',
                    choices=rad_choices,
                    help='choose radiation mechanism')

parser.add_argument('-emit',
                    action='store_true',
                    default=False,
                    help='enable photon emission')

parser.add_argument('-qed',
                    action='store_true',
                    default=False,
                    help='enable QED step')

parser.add_argument('-bwpp',
                    action='store_true',
                    default=False,
                    help='enable Breit-Wheeler pair production')

args = vars(parser.parse_args())

# Step 2. Set definitions and Makefile options based on above arguments

makefile_options = {}

if (args['user']):
  makefile_options['USER_FILE'] = args['user']
  makefile_options['USER_DIR'] = user_directory
else:
  makefile_options['USER_FILE'] = args['unit']
  makefile_options['USER_DIR'] = unit_directory

makefile_options['COMPILER_COMMAND'] = ''
makefile_options['COMPILER_FLAGS'] = ''
makefile_options['PREPROCESSOR_FLAGS'] = ''

# specific cluster:
specific_cluster = False
if args['perseus']:
  specific_cluster = True
  args['intel'] = True
  args['mpi08'] = True
  args['mpi'] = False
  args['ifport'] = True
  makefile_options['COMPILER_FLAGS'] += '-xCORE-AVX2 '

# compilation command
if args['hdf5']:
  makefile_options['COMPILER_COMMAND'] += 'h5pfc '
  makefile_options['PREPROCESSOR_FLAGS'] += '-DHDF5 '
else:
  if ((not args['mpi']) and (not args['mpi08'])):
    makefile_options['COMPILER_COMMAND'] += 'gfortran '
  else:
    makefile_options['COMPILER_COMMAND'] += 'mpif90 '
if args['ifport']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DIFPORT '

# mpi version
if args['mpi']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI '
elif args['mpi08']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI08 '

# debug
if args['debug'] and (not args['intel']):
  makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG -fcheck=all -fimplicit-none -fbacktrace '
if args['debug'] and args['intel']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG '
  makefile_options['COMPILER_FLAGS'] += '-traceback '

# compiler (+ vectorization etc)
if args['intel']:
  makefile_options['MODULE'] = '-module '
  makefile_options['COMPILER_FLAGS'] += '-O3 -DSoA -xHost -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto '
else:
  makefile_options['MODULE'] = '-J '
  makefile_options['COMPILER_FLAGS'] += '-O3 -DSoA -fwhole-program -mavx2 -fopt-info-vec -fopt-info-vec-missed -ftree-vectorizer-verbose=5 '

if args['1d']:
  makefile_options['EXE_NAME'] = 'tristan-mp1d'
  makefile_options['PREPROCESSOR_FLAGS'] += '-DoneD '
elif args['2d']:
  makefile_options['EXE_NAME'] = 'tristan-mp2d'
  makefile_options['PREPROCESSOR_FLAGS'] += '-DtwoD '
elif args['3d']:
  makefile_options['EXE_NAME'] = 'tristan-mp3d'
  makefile_options['PREPROCESSOR_FLAGS'] += '-DthreeD '

if args['absorb']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DABSORB '

# extra algorithms
if args['dwn']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DDOWNSAMPLING '
if args['alb'] and (not args['slb']):
  makefile_options['PREPROCESSOR_FLAGS'] += '-DALB '
if args['slb']:
  args['alb'] = False
  makefile_options['PREPROCESSOR_FLAGS'] += '-DSLB '
if args['gca']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DGCA '

# extra physics
if args['extfields']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DEXTERNALFIELDS '

if args['radiation'] != 'OFF':
  makefile_options['PREPROCESSOR_FLAGS'] += '-DRADIATION '

if 'sync' in args['radiation']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DSYNCHROTRON '
if 'ic' in args['radiation']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DINVERSECOMPTON '

if args['emit'] and (args['radiation'] != 'OFF'):
  makefile_options['PREPROCESSOR_FLAGS'] += '-DEMIT '

if args['qed']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DQED '

if args['bwpp']:
  makefile_options['PREPROCESSOR_FLAGS'] += '-DBWPAIRPRODUCTION '

makefile_options['PREPROCESSOR_FLAGS'] += '-DNGHOST=' + str(args['nghosts']) + ' '

# Step 3. Create new files, finish up
with open(makefile_input, 'r') as current_file:
  makefile_template = current_file.read()
for key, val in makefile_options.items():
  makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)
makefile_template = re.sub('# Template for ', '# ', makefile_template)
with open(makefile_output, 'w') as current_file:
  current_file.write(makefile_template)

# Finish with diagnostic output
print('==============================================================================')
print('Your TRISTAN distribution has now been configured with the following options:')
if (specific_cluster):
  if (args['perseus']):
    print('  Cluster configurations:  `Perseus`' )

print('SETUP ........................................................................')
print('  Userfile:                ' + makefile_options['USER_FILE'])
print('  Dim:                     ' + ('1D' if args['1d'] else ('2D' if args['2d'] else ('3D' if args['3d'] else 'None'))))
print('  # of ghost zones:        ' + str(args['nghosts']))
print('  Load balancing:          ' + ('adaptive' if args['alb'] else ('static' if args['slb'] else 'OFF')))
print('  Particle downsampling:   ' + ('ON' if args['dwn'] else 'OFF'))
print('  Particle pusher:         ' + ('Boris/GCA' if args['gca'] else 'Boris'))

print('PHYSICS ......................................................................')
print('  External fields:         ' + ('ON' if args['extfields'] else 'OFF'))
print('  Absorbing boundaries:    ' + ('ON' if args['absorb'] else 'OFF'))
print('  Cooling:                 ' + args['radiation'])
print('  Photon emission          ' + ('ON' if args['emit'] else 'OFF'))
print('  QED step                 ' + ('ON' if args['qed'] else 'OFF'))
print('  BW pair production       ' + ('ON' if args['bwpp'] else 'OFF'))

print('TECHNICAL ....................................................................')

print('  Compiler:                ' + ('intel' if args['intel'] else 'gcc'))
print('  Debug mode:              ' + ('ON' if args['debug'] else 'OFF'))
print('  Output:                  ' + ('HDF5' if args['hdf5'] else 'binary'))
print('  MPI version:             ' + ('old' if not args['mpi08'] else 'MPI_08'))
print('  `IFPORT` mkdir:          ' + ('ON' if args['ifport'] else 'OFF'))

print('==============================================================================')

print('  Compilation command:     ' + makefile_options['COMPILER_COMMAND'] \
    + makefile_options['PREPROCESSOR_FLAGS'] + makefile_options['COMPILER_FLAGS'])
