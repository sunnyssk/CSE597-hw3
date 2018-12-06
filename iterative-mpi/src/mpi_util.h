#ifndef _MPI_UTIL_H_
#define _MPI_UTIL_H_

#include <mpi.h>

// Synchronizes elements in an array with length of (size) distributed on different processes.
void MPIArraySync (double * array, int size, int mpi_size, int mpi_rank);

// Gathers the values in (target) on all processes and finds the maximum. Target will be modified to this maximum value.
void MPIGatherMax (double * target, int mpi_size, int mpi_rank);

// Starts MPI, rewrites mpi_size and mpi_rank.
void MPIStart (int & mpi_size, int & mpi_rank);

#endif /* _MPI_UTIL_H_ */