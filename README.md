# TRISTAN v2

#### Configuring & making
To configure simply type the following command
```bash
$ python configure.py [-FLAGS]
```
`[-FLAGS]` will be constantly updated. List of flags available now:

- `-3d`: enable 3D
- `-debug`: enable debug mode
- `-dprec`: use double precision (i.e., `real*8`)
- `-hdf5`: enable `hdf5` and compile with `h5fpc` (otherwise compiles with `gfortran`)
- `--user=[USER_FILE]`: name of the user file from `user/` directory (without the extension)
- _... more to come_

`Makefile` will be generated in the main directory (from `Makefile.in`). Now the code can be compiled with `make all` or cleaned with `make clean`.

#### Running
Executable is generated in `exec` directory named either `tristan-mp2d` or `tristan-mp3d`. Simply run with
```bash
$ exec/tristan-mp2d -i [input_file_name] -o [output_dir_name] -r [restart_dir_name]
```

#### Code structure
_TO BE ADDED_

#### Coding style advices
Please read this carefully.

1. We use JS style for naming variables and functions:
    - variables are named lowercase with underscores (`_`) if necessary, e.g.:
        - `my_new_var`;
    - function names start lowercase without underscores and can be continued uppercase, e.g.:
        - `myNewFunc()`;
    - modules always start with `m_` and are always lowercase without underscores, e.g.:
        - `m_mynewmodule`.
2. Remember, for fortran `THIS` and `tHIs` and `this` are the same.
3. Use indentations when entering loops, functions, conditional statements etc! Standard for this code is two spaces per one indent.
4. `#ifdef`-s and other precompiler statements can also be indented(!), since now `Makefile` precompiles the code with the new `cpp`.
5. Please use spaces and brackets to make code more readable:
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
