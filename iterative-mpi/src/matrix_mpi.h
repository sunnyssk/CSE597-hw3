#ifndef _MATRIX_MPI_H_
#define _MATRIX_MPI_H_

#include <iostream>
#include <cstring>
#include "mpi_util.h"

// Matrices are stored piecewise in processes, with (slice_rows_) rows each piece.
// Data is stored with column number changing in the most rapid manner.

template <typename T> class MMat {
public:
    MMat (int rows, int cols);
    ~MMat () { delete[] buffer_; }

    static void Init (int mpi_size, int mpi_rank);
    static void Clean ();

    int Rows() const { return rows_; }
    int Cols() const { return cols_; }
    int SliceRows() const { return slice_rows_; }
    int SliceOffset() const { return slice_offset_; }
    T * Buffer() { return buffer_; }

    // Gets the element in current process.
    T& Elem (int slice_row, int col) { return buffer_[slice_row * cols_ + col]; }
    const T& Elem (int slice_row, int col) const { return buffer_[slice_row * cols_ + col]; }
    // [Aligned] Gets the value of a certain element in matrix and broadcast it to all processes.
    T GetVal (int row, int col) const;
    
    // [Aligned] Print the matrix to the output stream specified.
    void Print (std::ostream& output_stream) const;

protected:
    static int mpi_size_;
    static int mpi_rank_;
    static MPI_Datatype mpi_type_;
    static bool initialized_;

    int rows_;
    int cols_;
    int slice_rows_;
    int slice_offset_;
    T* buffer_;
};

template class MMat<int>;
template class MMat<float>;
template class MMat<double>;

typedef MMat<int> MMatI;
typedef MMat<float> MMatF;
typedef MMat<double> MMatD;

#endif /* _MATRIX_MPI_H_ */