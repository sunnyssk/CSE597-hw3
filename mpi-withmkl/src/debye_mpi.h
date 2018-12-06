#ifndef _DEBYE_MPI_H_
#define _DEBYE_MPI_H_

#include <cmath>
#include <cstdio>
#include <cstring>
#include "mpi_util.h"
#include "matrix_mpi.h"
#include "mkl.h"

#ifndef MAX_ITER_NUM
    #define MAX_ITER_NUM (10000)
#endif

// Each single process has a unique copy of the field. To sync with the other processes,
// use MPIValueSync() directly. This will broadcast the certain part to all other processes.
class MField3D{
public:
    MField3D (int nx, int ny, int nz, double dx, double dy, double dz, int mpi_size, int mpi_rank);
    ~MField3D () { if(buffer_ != nullptr) delete[] buffer_; }

    int Nx () const { return nx_; }
    int Ny () const { return ny_; }
    int Nz () const { return nz_; }
    int N () const { return n_; }
    double Dx () const { return dx_; }
    double Dy () const { return dy_; }
    double Dz () const { return dz_; }

    double operator () (int x, int y, int z) const { return buffer_[(z * ny_ + y) * nx_ + x]; }
    double & operator () (int x, int y, int z) { return buffer_[(z * ny_ + y) * nx_ + x]; }

    // Applies Dirichlet's condition onto the field.
    void ApplyDirichletCond ();
    // [Aligned] Synchronizes values of the field distributed on different processes.
    void MPIValueSync ();
    // [Single process] Writes field into the file with specified filename.
    void WriteField (char const * filename);

protected:
    int nx_;
    int ny_;
    int nz_;
    int n_;
    int mpi_size_;
    int mpi_rank_;

    double dx_;
    double dy_;
    double dz_;
    double * buffer_;
};

class MDebyeSolver {
public:
    MDebyeSolver (int mpi_size, int mpi_rank) : mpi_size_(mpi_size), mpi_rank_(mpi_rank), pAmat_(nullptr), pfield_(nullptr) { }
    ~MDebyeSolver () { if(pAmat_ != nullptr) delete pAmat_; if(pfield_ != nullptr) delete[] pfield_; }

    // Generates the system matrix used in Jacobi iterative method.
    void GenerateSolverMatrix (MField3D const & rhs, double debye_length);
    // Carries out the Jacobi iterative solving and returns the iteration count. On output, the field copies are already the same on every process.
    int JacobiIterativeSolve (double err_threshold, MField3D & res_container, double * iter_err_array, int max_iter_num = MAX_ITER_NUM);
    // Inputs the right hand side vector for the linear problem.
    void RhsInput (MField3D const & field);
    
protected:
    int mpi_size_;
    int mpi_rank_;

    MMatD * pAmat_;
    double * pfield_;

    void MatMul(const double * x_prev, int M, int N, int slice_offset, double * x_next);
};

#endif /* _DEBYE_MPI_H_ */