#include "matrix_mpi.h"

template <typename T> int MMat<T>::mpi_rank_;
template <typename T> int MMat<T>::mpi_size_;
template <typename T> MPI_Datatype MMat<T>::mpi_type_;
template <typename T> bool MMat<T>::initialized_ = false;

template <typename T> MMat<T>::MMat (int rows, int cols) : 
rows_(rows), cols_(cols), slice_rows_((rows - 1) / mpi_size_ + 1), slice_offset_(slice_rows_ * mpi_rank_), buffer_(nullptr) {
    slice_rows_ = ((mpi_rank_ + 1) * slice_rows_ > rows) ? (rows - slice_rows_ * mpi_rank_) : (slice_rows_);
    slice_rows_ = slice_rows_ > 0 ? slice_rows_ : 0;
    buffer_ = new T[slice_rows_ * cols_];
}

template <typename T> void MMat<T>::Clean () {
    if (!initialized_) return;
    // else :
    MPI_Type_free(&mpi_type_);
    initialized_ = false;
}

template <typename T> T MMat<T>::GetVal (int row, int col) const {
    MPI_Barrier(MPI_COMM_WORLD);
    T retval;
    int slice_interval = (slice_offset_ == 0) ? slice_rows_ : (slice_offset_ / mpi_rank_), root = row / slice_interval;
    if (root == mpi_rank_) {
        retval = buffer_[(row - slice_offset_) * cols_ + col];
    }
    std::cout << "Root on # " << mpi_rank_ << ": " << root << std::endl;
    MPI_Bcast(&retval, 1, mpi_type_, root, MPI_COMM_WORLD);
    return retval;
}

template <typename T> void MMat<T>::Init (int mpi_size, int mpi_rank) {
    if (initialized_) return;
    // else :
    mpi_size_ = mpi_size;
    mpi_rank_ = mpi_rank;
    initialized_ = true;
    MPI_Type_contiguous(sizeof(T), MPI_CHAR, &mpi_type_);
    MPI_Type_commit(&mpi_type_);
}

template <typename T> void MMat<T>::Print (std::ostream& output_stream) const {
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank_ == 0) {
        T * res = new T[rows_ * cols_];
        memcpy(res, buffer_, sizeof(T) * slice_rows_ * cols_);
        int restrows = rows_;
        for (int i = 1; i < mpi_size_; i++) {
            restrows -= slice_rows_;
            int srcrows = (restrows >= slice_rows_) ? slice_rows_ : restrows;
            if(srcrows > 0) MPI_Recv(res + i * slice_rows_ * cols_, srcrows * cols_, mpi_type_, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        output_stream << "The matrix is:\n";
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) output_stream << res[i * cols_ + j] << " ";
            output_stream << "\n";
        }
        delete[] res;
    } else {
        if (slice_rows_ > 0) MPI_Send(buffer_, slice_rows_ * cols_, mpi_type_, 0, mpi_rank_, MPI_COMM_WORLD);
    }
}