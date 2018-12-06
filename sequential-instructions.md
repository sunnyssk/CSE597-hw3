## How to Compile on ACI

Switch to the root folder of the homework. Use the commands below:

    module load gcc/5.3.1
    make

Or you could use `make program` instead of `make`.

You should see "Make completed." if everything goes on well.

## Usage

Before executing the program, please make sure that you have a folder named "output" under the path you would like to call the program. If not, please create one.

Then you can run the program by specifying the input file:

    HW_ROOT_FOLDER/bin/LUDecomposition HW_ROOT_FOLDER/input/size.txt

The input file should contain only three integers (Nx, Ny, Nz), specifying the number of grids on each dimension.

## Expected Output

The program will show basic time & spatial consumptions on screen when running, along with the specific steps taken in solving the problems. The errors in iterative solver steps will be saved in `output/error-XXX.txt`, and results are given in `potential-XXX.txt` files. These output files could be processed and visualized with the iPython-Notebook files in the homework root folder.
