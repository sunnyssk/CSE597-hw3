#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mpi.h>
#include "mpi_util.h"
#include "matrix_mpi.h"
#include "debye_mpi.h"

int main (int argc, char** argv) {
    int mpi_size = 0, mpi_rank = 0;
    MPIStart(mpi_size, mpi_rank);
    MMatD::Init(mpi_size, mpi_rank);
    
    int nx = 20, ny = 20, nz = 20;
    double dl = 1E-9, debye_length = 2E-9;
    double ee = 1.60217662E-19, e0 = 8.854187817E-12;
    MField3D rhs(nx, ny, nz, dl, dl, dl, mpi_size, mpi_rank), potential(nx, ny, nz, dl, dl, dl, mpi_size, mpi_rank);
    MDebyeSolver msolve(mpi_size, mpi_rank);
    double *error_array = new double[MAX_ITER_NUM];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                potential(i, j, k) = 0.0;
                rhs(i, j, k) = 0.0;
            }
        }
    }

    rhs(nx / 2, ny / 2, nz / 2) = -10 * ee / (e0 * dl * dl * dl);               // -charge_density / epsilon
    msolve.GenerateSolverMatrix(rhs, debye_length);
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_rank == 0) std::cout << "Solver matrix generated." << std::endl;
    msolve.RhsInput(rhs);
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_rank == 0) std::cout << "RHS input completed." << std::endl;
    msolve.JacobiIterativeSolve(1E-5, potential, error_array);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) potential.WriteField("output/potential.txt");

    delete[] error_array;
    MMatD::Clean();
    MPI_Finalize();
    return 0;
}