# Readme: HW03

## Author Info

Author     : Yueze Tan

Email      : yut75@psu.edu

Last update: 2018/12/05

License    : GNU (See `license.txt` for details. A copy of GPL is provided as in `gpl.txt`.)

## How to Compile on ACI

Depending on whether you want to use the profiling tool or not, the compiling procedures are different.

### Compile without profile

If you want to use the program without profiling, please switch to the folder `mpi-withmkl`, then open `Makefile`, uncomment the line with `CC=mpicxx`, and comment the line `CC=tau_cxx.sh`.

Now save and quit the editor. Load the modules `gcc/5.3.1`, `mkl` and `openmpi` (with the default version) and modify the `MKL_INCLUDE_PATH` and `MKL_LIB_PATH` in Makefile to the expected paths. Then execute

    make clean
    make

Or you could use `make program` instead of `make`.

You should see "Make completed." if everything goes on well.

### Compile with profile

We use `tau` for profile. For the code is based on MPI, it is possible to use the profiling tool provided by the lecturer.

If you want to use the program with profiling, the you should uncomment the line `CC=tau_cxx.sh` in the Makefile, and comment `CC=mpicxx`. Except for loading `gcc/5.3.1`, `mkl` and `openmpi`, to compile the parallel code, please use

    module use /storage/work/a/awl5173/toShare/tauPdt/tau
    module load adamsTau_2.27

or any other version of `tau` with supporting of openmpi if you have one. The requests for MKL library and include paths are the same as compiling without profiling tool.

To start compiling, type

    make clean
    make

It is possible that the profiler generates much info and requires you to type ENTER for multiple times during the compilation. If you do not see any error reported by the compiler (rather than profiler), and the link appears to be successful, then the program will be ready to use.

## Usage

Before executing the program, please make sure that you have a folder named "output" under the path you would like to call the program. If not, please create one. Or the program will not be able to output the final results.

Then you can run the program with

`mpirun -np N HW_ROOT_FOLDER/mpi/bin/MKLJacobi`

with `N` being the number of processes you desire.

After execution, the profile will be generated if you have compiled with tau. In that case, you can see the profile with `pprof` or `paraprof` with the tau module.

## Expected Output

The program will show the specific iteration steps taken in solving the problems. The results are given in `potential.txt`. The output files could be processed and visualized with the iPython-Notebook file in the homework root folder.

## Instructions for Code without MKL

To compile and run the iterative solver codes without MKL support in `iterative-mpi` and `iterative-sequential`, see `iterative-instructions.md`.

## Instructions for Sequential Solvers

To compile and run the sequential codes for direct and iterative solvers in `direct-and-iterative`, see `sequential-instructions.md`.

## Write-up Compilation

The write-up is generated with [overleaf](https://v2.overleaf.com). A copy is already included in the folder `writeup`.

To generate a new copy, please upload the .tex and .bib file to the website, along with the three figures in the `output` folder.

Create an `output` and a `writeup` folder under one project, place all the figures under `output`, and place the other files under `writeup`.

Set main file of project as `writeup/main.tex`, then click the compile button to get the output.

## Acknowledgements

Please see the LaTeX write-up and its PDF output for details.

