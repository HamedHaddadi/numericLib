# ifndef JACOBI_H_
# define JACOBI_H_

# include "../utils/header.h"
# include "../arrays/matrix1D.h"
# include "../arrays/matrix2D.h"
# include "../utils/concepts.h"


template <typename T>
requires std::is_floating_point<T>::value && TwoD<Matrix2D<T>> && OneD<Matrix1D<T>>
Matrix1D<T> jacobi(const Matrix2D<T> & A, const Matrix1D<T> & b, int maxIter, double tol = 1e-5) {
    Matrix1D<T> xZero(b.getShape());
    Matrix1D<T> x(b.getShape());
    std::pair<std::size_t, std::size_t> matShape = A.getShape();
    T colSum_{};
    T norm_{};
    for (int iter = 1; iter < maxIter; ++iter) {
        for (std::size_t row = 0; row < matShape.first; ++row) {
            colSum_ = T(0);
            for (std::size_t col = 0; col < matShape.second; ++col) {
                if (col != row)
                    colSum_ += -1.0*A(row, col)*xZero(col);
            }
            x(row) = (1./A(row, row))*(colSum_ + b(row));
        }
        norm_ = euclidean(x, xZero);
        if (norm_ < tol) 
            break;
        xZero = x;
    }
    return x;
} 

template <typename U, typename W>
requires TwoD<Matrix2D<U>> && OneD<Matrix1D<W>>
Matrix1D<std::common_type_t<U, W> > jacobi(const Matrix2D<U> & A, const Matrix1D<W> & b, int maxIter, double tol = 1e-5) {
    typedef std::common_type_t<U, W> T; 
    Matrix1D<T> x(b.getShape());
    Matrix1D<T> xZero(b.getShape());
    std::pair<std::size_t, std::size_t> matShape = A.getShape();
    T colSum_{};
    T norm_{};
    for (int iter = 1; iter < maxIter; ++iter) {
        for (std::size_t row = 0; row < matShape.first; ++row) {
            colSum_ = T(0);
            for (std::size_t col = 0; col < matShape.second; ++col) {
                if (col != row)
                    colSum_ += -1.0*static_cast<T>(A(row, col))*xZero(col);
            }
            x(row) = (1./static_cast<T>(A(row, row)))*(colSum_ + b(row));
        }
        norm_ = euclidean(x, xZero);
        if (norm_ < tol) 
            break;
        xZero = x;
    }
    return x;    
} 



# endif 