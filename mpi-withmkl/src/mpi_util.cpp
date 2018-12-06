#include "mpi_util.h"

void MPIArraySync (double * array, int size, int mpi_size, int mpi_rank) {
    int mpi_interval = (size - 1) / mpi_size + 1, slice_rows = 0, offset = 0;
    for (int i = 0; i < mpi_size; i++) {
        slice_rows = ((i + 1) * mpi_interval > size) ? (size - mpi_interval * i) : (mpi_interval);
        slice_rows = slice_rows > 0 ? slice_rows : 0;
        if(slice_rows > 0){
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(array + offset, slice_rows, MPI_DOUBLE, i, MPI_COMM_WORLD);
        }
        offset += mpi_interval;
    }
}

void MPIGatherMax (double * target, int mpi_size, int mpi_rank) {
    double buffer = 0.0;
    for (int i = 0; i < mpi_size; i++) {
        if(i == mpi_rank) buffer = *target;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&buffer, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
        if(buffer > *target) *target = buffer;
    }
}

void MPIStart (int & mpi_size, int & mpi_rank) {
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
}