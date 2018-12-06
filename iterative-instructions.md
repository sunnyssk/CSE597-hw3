## How to Compile on ACI

Depending on whether you want to use the profiling tool or not, the compiling procedures are different.

### Compile without profile

If you want to use the program without profiling, then switch to the folder `sequential` or `mpi` (serial or parallel, depending on which one you want to use), open `Makefile`, uncomment the line with `CC=g++` or `CC=mpicxx`, and comment the line `CC=tau_cxx.sh`.

Now save and quit the editor. Load the modules `gcc/5.3.1` and `openmpi` (with the default version) and then type

    make clean
    make

Or you could use `make program` instead of `make`.

You should see "Make completed." if everything goes on well.

### Compile with profile

We use `tau` for profile. For the code is based on MPI, it is possible to use the profiling tool provided by the lecturer.

If you want to use the program with profiling, the you should uncomment the line `CC=tau_cxx.sh` in the Makefile, and comment `CC=g++` or `CC=mpicxx`. Except for loading `gcc/5.3.1` and `openmpi`, to compile the serial code, you will be needing `tau/2.2.7` on ACI; to compile the parallel code, please use

    module use /storage/work/a/awl5173/toShare/tauPdt/tau
    module load adamsTau_2.27

or any other version of `tau` with supporting of openmpi if you have one.

To start compiling, type

    make clean
    make

It is possible that the profiler generates much info and requires you to type ENTER for multiple times during the compilation. If you do not see any error reported by the compiler (rather than profiler), and the link appears to be successful, then the program will be ready for use.

## Usage

Before executing the program, please make sure that you have a folder named "output" under the path you would like to call the program. If not, please create one. Or the program will not be able to output the final results.

Then you can run the program with

`HW_ROOT_FOLDER/sequential/bin/MPIJacobi` (serial version)

or

`mpirun -np N HW_ROOT_FOLDER/mpi/bin/MPIJacobi` (parallel version)

with `N` being the number of processes you desire.

After execution, the profile will be generated if you have compiled with tau. In that case, you can see the profile with `pprof` or `paraprof` with the tau module.

## Expected Output

The program will show the specific iteration steps taken in solving the problems. The results are given in `potential.txt`. The output files could be processed and visualized with the iPython-Notebook file in the homework root folder.
