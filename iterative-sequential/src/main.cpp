#include <iostream>
#include <cstdlib>
#include <cmath>
#include <exception>
#include <stdio.h>

#include "matrix.h"
#include "debye.h"

int main (int argc, char** argv) {
    srand(time(NULL));
    int nx = 20, ny = 20, nz = 20;
    int n = nx * ny * nz;
    double ee = 1.60217662E-19, e0 = 8.854187817E-12;
    double debye_length = 2E-9;
    double * error_array = new double[MAX_ITER_NUM];
    Field3D rhs(nx, ny, nz, 1E-9), potential1(nx, ny, nz, 1E-9), potential2(nx, ny, nz, 1E-9), potential3(nx, ny, nz, 1E-9);
    for (int i = 0; i < 30; i++) rhs(rand() % (nx - 2) + 1, rand() % (ny - 2) + 1, rand() % (nz - 2) + 1) = -10 * ee / e0 * 1E27;       // -rho / epsilon_0

    DebyeSolver dsolve;
    dsolve.GenerateSolverMatrix(rhs, debye_length);
    dsolve.RhsInput(rhs);
    
    dsolve.JacobiIterativeSolve(1E-5, potential2, error_array);
    FILE *fout = fopen("output/potential.txt", "w");
    potential2.WriteField(fout);
    fclose(fout);

    delete[] error_array;
    return 0;
}
