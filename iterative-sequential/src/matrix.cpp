#include "matrix.h"

template <typename T> void Matrix<T>::Add (Matrix<T> const & opr2) {
    if (rows_ != opr2.rows_ || cols_ != opr2.cols_) throw "ERROR: Arithmatic operation on matrices with different sizes.";
    // else :
    int size = rows_ * cols_;
    for (int i = 0; i < size; i++) data_[i] += opr2.data_[i];
}

template <typename T> void Matrix<T>::Add (Matrix<T> const & opr1, Matrix<T> const & opr2, Matrix<T> & res) {
    if (opr1.rows_ != opr2.rows_ || opr1.cols_ != opr2.cols_) throw "ERROR: Arithmatic operation on matrices with different sizes.";
    else if (opr1.rows_ != res.rows_ || opr2.cols_ != res.cols_) throw "ERROR: Size of result container does not match the size of input.";
    // else:
    int size = opr1.rows_ * opr1.cols_;
    for (int i = 0; i < size; i++) res.data_[i] = opr1.data_[i] + opr2.data_[i];
}

template <typename T> void Matrix<T>::Copy (Matrix<T> const & src) {
    int size = src.rows_ * src.cols_;
    if (src.rows_ != rows_ || src.cols_ != cols_) {
        delete[] data_;                     // releases current data array first
        data_ = new T[size];
        rows_ = src.rows_;
        cols_ = src.cols_;
        mapping_only_ = false;
    } // else: reuse the existing memory
    for (int i = 0; i < size; i++) data_[i] = src.data_[i]; 
}

template <typename T> void Matrix<T>::Eye (int rows) {
    int size = rows * rows;
    if (rows != rows_ || rows != cols_) {
        delete[] data_;                     // releases current data array first
        data_ = new T[size];
        rows_ = rows;
        cols_ = rows;
        mapping_only_ = false;
    } else for (int i = 0; i < size; i++) this->data_[i] = T(0);
    for (int i = 0; i < rows; i++) (*this)(i, i) = T(1);
}

template <typename T> void Matrix<T>::LInverse (Matrix<T> & res_container) const {
    res_container.Copy(*this);
    for (int k = 0; k < rows_; k++) {
        T ratio = T(1) / res_container(k, k);
        for (int j = 0; j < k; j++) res_container(k, j) *= ratio;
        for (int i = k + 1; i < rows_; i++) {
            for (int j = 0; j < k; j++) res_container(i, j) -= res_container(i, k) * res_container(k, j);
            res_container(i, k) = -res_container(i, k) * ratio;
        }
        res_container(k, k) = ratio;
    }
}

template <typename T> void Matrix<T>::LSolve (Matrix<T> const & src, Matrix<T> & res_container) const {
    int n = this->Rows();
    for (int i = 0; i < n; i++) {
        T sum = T(0);
        for (int j = 0; j < i; j++) sum += (*this)(i, j) * res_container(j, 0);
        res_container(i, 0) = (src(i, 0) - sum) / (*this)(i, i);
    }
}

template <typename T> void Matrix<T>::LUDecomposition (Matrix<T> & l_container, Matrix<T> & u_container) const {
    if (rows_ != cols_) throw "ERROR: LU decomposition not applying on a square matrix.";
    // else:
    u_container.Copy(*this);
    l_container.Zeros(rows_, cols_);
    for (int k = 0; k < rows_; k++) {
        l_container(k, k) = T(1);
        for (int i = k + 1; i < rows_; i++) {
            T ratio = -u_container(i, k) / u_container(k, k);
            u_container(i, k) = T(0);
            l_container(i, k) = -ratio;
            for (int j = k + 1; j < cols_; j++) u_container(i, j) += ratio * u_container(k, j);
        }
    }
}

template <typename T> void Matrix<T>::Multiply (Matrix<T> const & opr1, Matrix<T> const & opr2, Matrix<T> & res) {
    if (opr1.cols_ != opr2.rows_ || opr1.rows_ != res.rows_ || opr2.cols_ != res.cols_) throw "ERROR: Matrix operand dimensions do not match.";
    // else:
    memset(res.data_, T(0), sizeof(T) * res.Size());
    for (int i = 0; i < opr1.rows_; i++)
        for (int j = 0; j < opr2.cols_; j++)
            for (int k = 0; k < opr1.cols_; k++) res(i, j) += opr1(i, k) * opr2(k, j);
}

template <typename T> void Matrix<T>::Multiply (T const & opr) {
    int size = rows_ * cols_;
    for (int i = 0; i < size; i++) data_[i] *= opr;
}

template <typename T> void Matrix<T>::PLUDecomposition (Matrix<T> & p_container, Matrix<T> & l_container, Matrix<T> & u_container) const {
    if (rows_ != cols_) throw "ERROR: PLU decomposition not applying on a square matrix.";
    // else:
    Matrix<T> u_permutate(rows_, cols_), l_permutate(rows_, cols_);
    u_permutate.Copy(*this);
    l_permutate.Zeros(rows_, cols_);
    p_container.Zeros(rows_, cols_);
    int* p = new int[rows_];                   // array records row permutations
    for (int i = 0; i < rows_; i++) p[i] = i;
    for (int k = 0; k < rows_; k++) {
        int pivot = k, tmp = 0;
        for (int i = k; i < rows_; i++)                         // find the pivot
            if (fabs(u_permutate(p[pivot], k)) < fabs(u_permutate(p[i], k))) pivot = i;
        tmp = p[k];
        p[k] = p[pivot];
        p[pivot] = tmp;

        l_permutate(p[k], k) = T(1);
        for (int i = k + 1; i < rows_; i++) {
            T ratio = -u_permutate(p[i], k) / u_permutate(p[k], k);
            u_permutate(p[i], k) = T(0);
            l_permutate(p[i], k) = -ratio;
            for (int j = k + 1; j < cols_; j++) u_permutate(p[i], j) += ratio * u_permutate(p[k], j);
        }
    }
    for (int i = 0; i < rows_; i++) {
        p_container(i, p[i]) = T(1);
        for (int j = 0; j < cols_; j++) {
            l_container(i, j) = l_permutate(p[i], j);
            u_container(i, j) = u_permutate(p[i], j);
        }
    }
    delete[] p;
}

template <typename T> void Matrix<T>::Print (std::ostream & output_stream) const {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_ - 1; j++) {
            output_stream << data_[i + j * rows_] << ", ";
        }
        output_stream << data_[i + (cols_ - 1) * rows_] << "\n";
    }
}

template <typename T> void Matrix<T>::Rand (T const & lower_bound, T const & upper_bound) {
    for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++) (*this)(i, j) = (T(rand()) * (upper_bound - lower_bound)) / RAND_MAX + lower_bound;
}

template <typename T> void Matrix<T>::Reshape (int new_rows, int new_cols) {
    if (new_rows * new_cols != rows_ * cols_) throw "ERROR: Matrix size changed during reshape.";
    // else:
    rows_ = new_rows;
    cols_ = new_cols;
}

template <typename T> void Matrix<T>::Subtract (Matrix<T> const & opr2) {
    if (rows_ != opr2.rows_ || cols_ != opr2.cols_) throw "ERROR: Arithmatic operation on matrices with different sizes.";
    // else :
    int size = rows_ * cols_;
    for (int i = 0; i < size; i++) data_[i] -= opr2.data_[i];
}

template <typename T> void Matrix<T>::Subtract (Matrix<T> const & opr1, Matrix<T> const & opr2, Matrix<T> & res) {
    if (opr1.rows_ != opr2.rows_ || opr1.cols_ != opr2.cols_) throw "ERROR: Arithmatic operation on matrices with different sizes.";
    else if (opr1.rows_ != res.rows_ || opr2.cols_ != res.cols_) throw "ERROR: Size of result container does not match the size of input.";
    // else:
    int size = opr1.rows_ * opr1.cols_;
    for (int i = 0; i < size; i++) res.data_[i] = opr1.data_[i] - opr2.data_[i];
}

template <typename T> void Matrix<T>::SwapRows (int row1, int row2) {
    T tmp;
    for (int i = 0; i < cols_; i++) {
        int ind1 = i * rows_ + row1, ind2 = i * rows_ + row2;
        tmp = data_[ind1];
        data_[ind1] = data_[ind2];
        data_[ind2] = tmp;
    }
}

template <typename T> void Matrix<T>::SwapCols (int col1, int col2) {
    T tmp;
    for (int i = 0; i < rows_; i++) {
        int ind1 = col1 * rows_ + i, ind2 = col2 * rows_ + i;
        tmp = data_[ind1];
        data_[ind1] = data_[ind2];
        data_[ind2] = tmp;
    }
}

template <typename T> Matrix<T> Matrix<T>::Transpose() {
    Matrix<T> result(cols_, rows_);
    for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++) result(j, i) = (*this)(i, j);
    return result;
}

template <typename T> void Matrix<T>::UInverse (Matrix<T> & res_container) const {
    res_container.Copy(*this);
    for (int k = rows_ - 1; k >= 0; k--) {
        T ratio = T(1) / res_container(k, k);
        for (int j = k + 1; j < cols_; j++) res_container(k, j) *= ratio;
        for (int i = k - 1; i >= 0; i--) {
            for (int j = k + 1; j < cols_; j++) res_container(i, j) -= res_container(i, k) * res_container(k, j);
            res_container(i, k) = -res_container(i, k) * ratio;
        }
        res_container(k, k) = ratio;
    }
}

template <typename T> void Matrix<T>::USolve (Matrix<T> const & src, Matrix<T> & res_container) const {
    int n = this->Rows();
    for (int i = n - 1; i >= 0; i--) {
        T sum = T(0);
        for (int j = i + 1; j < n; j++) sum += (*this)(i, j) * res_container(j, 0);
        res_container(i, 0) = (src(i, 0) - sum) / (*this)(i, i);
    }
}

template <typename T> void Matrix<T>::Zeros (int rows, int cols) {
    int new_size = rows * cols;
    if (rows != rows_ || cols != cols_) {
        delete[] data_;
        rows_ = rows;
        cols_ = cols;
        mapping_only_ = false;
        data_ = new T[new_size];
    } else for (int i = 0; i < new_size; i++) data_[i] = T(0);
}
