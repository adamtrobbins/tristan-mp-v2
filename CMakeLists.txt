cmake_minimum_required(VERSION 3.16)

# hdf5 enabled by default [hardcoded for now]
set(CMAKE_Fortran_COMPILER $ENV{HDF5_ROOT}/bin/h5pfc)

project(tristan-v2 VERSION 2.8
        DESCRIPTION "Relativistic cartesian particle-in-cell code"
        LANGUAGES Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(HDF5 REQUIRED COMPONENTS Fortran)

set(DEBUG "OFF" CACHE BOOL "Debug mode")

set(DEFINITIONS "\
  -DSTR_MAX=280\
  -DTINYXYZ=1e-6\
  -DTINYREAL=1e-3\
  -DTINYFLD=1e-8\
  -DTINYWEI=1e-6\
  -DM_PI=3.141592653589793\
  -DVEC_LEN=16")

  # !TODO: debug leves & arch specific optimizations
if (${DEBUG} STREQUAL "ON")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDEBUG")
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -fimplicit-none -fbacktrace")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback")
  endif()
else()
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -DSoA -fwhole-program -mavx2 -fopt-info-vec -fopt-info-vec-missed -ftree-vectorizer-verbose=5 -ffree-line-length-512 ")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -DSoA -xCORE-AVX512 -ipo -qopenmp-simd -qopt-report=5 -qopt-zmm-usage=high -qopt-streaming-stores auto")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 10397 -diag-disable 10346 -diag-disable 8100 -diag-disable 6178")
  endif()
endif()

# -- options --
# architecture
option(MPI08 "Use MPI08" ON)
option(MPINONBLOCK "Use non-blocking MPI" OFF)
option(LOWMEM "Use low memory mode" OFF)
option(IFPORT "Use ifport" OFF)
# optimizations
option(DOWNSAMPLING "Use particle downsampling" OFF)
option(SLB "Use static load balancing" OFF)
option(ALB "Use adaptive load balancing" OFF)
option(PRTLPAYLOADS "Use particle payloads" OFF)
# output
option(SERIALOUTPUT "Output in serial" OFF)
option(USROUTPUT "Enable custom user data outputting" OFF)
# physics
option(EXTERNALFIELDS "Use external fields" OFF)
option(ABSORB "Use absorbing boundary conditions" OFF)
option(VAY "Use Vay current deposition" OFF)
option(BLINNE "Use Blinne field solvers" OFF)
option(EMIT "Enable photon emission" OFF)
option(QED "Enable QED" OFF)
option(COMPTONSCATTERING "Enable Compton scattering" OFF)
option(BWPAIRPRODUCTION "Enable Breit-Wheeler pair production" OFF)
option(PAIRANNIHILATION "Enable pair annihilation" OFF)

set(PRECISION "single" CACHE STRING "Precision of the code")
set(NGHOST "2" CACHE STRING "Number of ghost cells")
set(RADIATION "OFF" CACHE STRING "Radiation drag")

# -- additional configs for each option --
# !TODO
add_definitions(-D NGHOST=${NGHOST})

add_definitions(${DEFINITIONS} -D MPI08)
add_subdirectory(src)
