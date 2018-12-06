# Readme: HW02

## Author Info

Author     : Yueze Tan

Email      : yut75@psu.edu

Last update: 2018/10/26

License    : GNU (See `license.txt` for details. A copy of GPL is provided as in `gpl.txt`.)

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

## Write-up Compilation

The write-up is generated with [overleaf](https://v2.overleaf.com). A copy is already included in the folder `writeup`.

To generate a new copy, please upload the .tex and .bib file to the website, along with the three figures in the `output` folder.

Create an `output` and a `writeup` folder under one project, place all the figures under `output`, and place the other files under `writeup`.

Set main file of project as `writeup/main.tex`, then click the compile button to get the output.

## Acknowledgements

Please see the LaTeX write-up and its PDF output for details.

