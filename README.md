# TRISTAN v2

#### Configuring & making
To configure simply type the following command
```bash
$ python configure.py [-FLAGS]
```
Supported `[-FLAGS]` will be constantly updated. List of flags available now:

- `-3d`: enable 3D;
- `-intel`: compile with intel compatibility;
- `-debug`: enable debug mode (enables custom `-DDEBUG` macros flag, `traceback`, `qopt` reports for intel compilers etc);
- `--nghosts=[NUM_GHOSTS]`: specify the number of ghost zones, will compile with `-DNGHOST=[NUM_GHOSTS]`;
- `--user=[USER_FILE]`: name of the user file from `user/` directory (without the extension);
- `-hdf5`: enable `hdf5` and compile with `h5fpc`;
- `-ifport`: handles `mkdir` commands, some systems do not support this;
- _... more to come_.


`Makefile` will be generated in the main directory (from `Makefile.in`). Now the code can be compiled with `make all` or cleaned with `make clean`.

#### Prerequisites
For MPI we use `mpi_f08` standard.

#### Running
Executable is generated in `exec` directory named either `tristan-mp2d` or `tristan-mp3d`. Simply run with
```bash
$ exec/tristan-mp2d -i [input_file_name] -o [output_dir_name] -r [restart_dir_name]
```
or if compiled with `mpif90` one can run in MPI
```bash
$ mpiexec -np [NPROC] exec/tristan-mp2d -i [input_file_name] -o [output_dir_name] -r [restart_dir_name]
```

#### Done so far
- Logistics: `configure.py`, read input, initialize everything, etc;
- CPU & "meshblocks" handling;
- particle tiles (size is configurable from the input);
- species + particles: array of structures (species) of arrays (`xi`, `yi`, `zi`, etc);
- particle pusher (aligned) & particle exchange;
- fields & field exchange (ghost zones);
- binary/hdf5 output of particles, fields, spectra + `python` library to read and create a readable dictionary (and convert to `hdf5` if necessary).

#### ToDo
- restart files;
- "history" file;
- absorbing (radiation) boundaries;
- expanding boundaries;
- static/adaptive load balancing;
- pair-production/annihilation/IC routine;
- make it work in 1D;
- generalize for non-MPI.

#### Code structure
_TO BE ADDED_

#### Coding style advices
Please read this carefully.

1. We use the following style for naming variables and functions:
    - variables are named lowercase with underscores (`_`) if necessary, e.g.:
        - `my_new_var`;
    - function names start lowercase without underscores and can be continued uppercase, e.g.:
        - `myNewFunc()`;
    - modules always start with `m_` and are always lowercase without underscores, e.g.:
        - `m_mynewmodule`.
2. Remember, for fortran `THIS` and `tHIs` and `this` are the same.
3. Use indentations when entering loops, functions, conditional statements etc! Standard for this code is two spaces per one indent.
4. `#ifdef`-s and other precompiler macros can also be indented(!), as now `Makefile` precompiles the code with the new `cpp` compiler before compiling with fortran.
5. Please use spaces and brackets to make the code more readable:
    - in arithmetic statements, e.g.:
        - `2 + 3 * 5 / (2 + a)` instead of `2+3*5/(2+a)`;
    - in function arguments, e.g.:
        - `myFunction(2, var, "my string")` instead of `myFunction(2,var,"my string")`;
    - in assignments, e.g.:
        - `my_var = my_other_var` instead of `my_var=my_other_var`;
    - in conditional statements, e.g.:
        - `if ((one .eq. 1) .or. (two .eq. 2))` instead of `if(one.eq.1 .or. two.eq.2)`.
6. Try to use standard Fortran `.eq.`, `.lt.` etc instead of `==`, `<` etc.
7. Use `interface`-s and `module procedure` to overload some of the auxiliary functions and make them accept generic arguments.
8. Do __not__ use `goto` statements.

All these advices and rules are simply to make the code look visually comprehensible and help those working with it enjoy the time. Thank you for following these advices.
