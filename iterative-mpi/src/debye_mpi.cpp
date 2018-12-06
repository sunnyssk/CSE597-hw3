#include "debye_mpi.h"

MField3D::MField3D (int nx, int ny, int nz, double dx, double dy, double dz, int mpi_size, int mpi_rank) :
nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), mpi_size_(mpi_size), mpi_rank_(mpi_rank) {
    n_ = nx_ * ny_ * nz_;
    buffer_ = new double [n_];
}

void MField3D::ApplyDirichletCond () {
    int nx = nx_, ny = ny_, nz = nz_;
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            (*this)(0, j, k) = 0.0;
            (*this)(nx - 1, j, k) = 0.0;
        }
    }
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            (*this)(i, 0, k) = 0.0;
            (*this)(i, ny - 1, k) = 0.0;
        }
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            (*this)(i, j, 0) = 0.0;
            (*this)(i, j, nz - 1) = 0.0;
        }
    }
}

void MField3D::MPIValueSync () {
    MPIArraySync(buffer_, n_, mpi_size_, mpi_rank_);
}

void MField3D::WriteField (char const * filename) {
    FILE *fout = fopen(filename, "w");
    if (!fout) std::cout << "Cannot create the output file, please make sure you have access to writing at the specified path." << std::endl;
    else {
        fprintf(fout, "%d %d %d\n", nx_, ny_, nz_);
        for (int i = 0; i < nx_; i++)
            for (int j = 0; j < ny_; j++)
                for (int k = 0; k < nz_; k++) fprintf(fout, "%lf\n", (*this)(i, j, k));
    }
}

void MDebyeSolver::GenerateSolverMatrix (MField3D const & rhs, double debye_length) {
    int nx = rhs.Nx(), ny = rhs.Ny(), nz = rhs.Nz();
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nltrim * nztrim;
    double dx = rhs.Dx(), dy = rhs.Dy(), dz = rhs.Dz();
    double dx2 = dx * dx, dy2 = dy * dy, dz2 = dz * dz, ldi2 = 1 / (debye_length * debye_length);
    if (pAmat_) delete pAmat_;
    pAmat_ = new MMatD(ntrim, ntrim);
    MMatD &A = *pAmat_;
    int row_offset = A.SliceOffset(), row_terminate = row_offset + A.SliceRows();

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int k = 1; k <= nztrim; k++) {
        if (TRI(0, 0, k + 1) <= row_offset || TRI(0, 0, k) > row_terminate) continue;
        for (int j = 1; j <= nytrim; j++) {
            if (TRI(0, j + 1, k) <= row_offset || TRI(0, j, k) > row_terminate) continue;
            for (int i = 1; i <= nxtrim; i++) {
                if (TRI(i + 1, j, k) <= row_offset || TRI(i, j, k) > row_terminate) continue;
                A.Elem(TRI(i, j, k) - row_offset, TRI(i, j, k)) = -(2 / dx2 + 2 / dy2 + 2 / dz2 + ldi2);
                if(i > 1) A.Elem(TRI(i, j, k) - row_offset, TRI(i - 1, j, k)) = 1 / dx2;
                if(j > 1) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j - 1, k)) = 1 / dy2;
                if(k > 1) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j, k - 1)) = 1 / dz2;
                if(i < nxtrim) A.Elem(TRI(i, j, k) - row_offset, TRI(i + 1, j, k)) = 1 / dx2;
                if(j < nytrim) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j + 1, k)) = 1 / dy2;
                if(k < nztrim) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j, k + 1)) = 1 / dz2;
            }
        }
    }
    #undef TRI
}

int MDebyeSolver::JacobiIterativeSolve (double err_threshold, MField3D & res_container, double * iter_err_array, int max_iter_num) {
    int nx = res_container.Nx(), ny = res_container.Ny(), nz = res_container.Nz();
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    int iter_cnt = 0;
    int slice_rows = pAmat_->SliceRows(), slice_offset = pAmat_->SliceOffset();

    double *aux_vector1 = new double[ntrim], *aux_vector2 = new double[ntrim];
    double *x_prev = aux_vector1, *x_next = aux_vector2;
    double errmax = 0.0;

    // Assign the initial guess to the input vector
    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) x_prev[TRI(i, j, k)] = res_container(i, j, k);
    #undef TRI

    MPI_Barrier(MPI_COMM_WORLD);
    // Iterative solver
    do {
        errmax = 0.0;
        if (iter_cnt++ >= max_iter_num - 1) {
            std::cout << "Iteration number exceeds limit! Exit automatically." << std::endl;
            break;
        }
        for (int i = slice_offset; i < slice_offset + slice_rows; i++) x_next[i] = pfield_[i];

        MatMul(x_prev, slice_rows, ntrim, slice_offset, x_next);
        for (int i = slice_offset; i < slice_offset + slice_rows; i++) {
            int slice_i = i - slice_offset;
            x_next[i] += pAmat_->Elem(slice_i, i) * x_prev[i];
            x_next[i] /= pAmat_->Elem(slice_i, i);
            double newerr = fabs(x_next[i] - x_prev[i]);
            errmax = errmax > newerr ? errmax : newerr;
        }
        
        MPIArraySync(x_next, ntrim, mpi_size_, mpi_rank_);
        MPIGatherMax(&errmax, mpi_size_, mpi_rank_);
        double *ptmp = x_prev;                                                  // swap two vectors
        x_prev = x_next;
        x_next = ptmp;
        if(mpi_rank_ == 0) std::cout << "Iteration round #" << iter_cnt << ", maximum error: " << errmax << "\n";
        iter_err_array[iter_cnt - 1] = errmax;
    } while (errmax >= err_threshold || aux_vector1 == x_next);                 // only exits iteration when the original field is updated

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) res_container(i, j, k) = x_next[TRI(i, j, k)];
    #undef TRI

    // char *buffer = new char[256];
    // MemSizeOutput(buffer);
    // delete[] buffer;
    delete[] aux_vector1;
    delete[] aux_vector2;
    return iter_cnt;
}

void MDebyeSolver::MatMul (const double * x_prev, int M, int N, int slice_offset, double * x_next) {
    for (int i = slice_offset; i < slice_offset + M; i++)
        for (int j = 0; j < N; j++) x_next[i] -= pAmat_->Elem(i - slice_offset, j) * x_prev[j];
}

void MDebyeSolver::RhsInput (MField3D const & field) {
    int nx = field.Nx(), ny = field.Ny(), nz = field.Nz(), n = nx * ny * nz;
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    if (pfield_) delete[] pfield_;
    pfield_ = new double[ntrim];

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) pfield_[TRI(i, j, k)] = field(i, j, k);
    #undef TRI
}
