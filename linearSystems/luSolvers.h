# ifndef LUDECOMPOSITION_H_
# define LUDECOMPOSITION_H_

# include "../utils/header.h"
# include "../utils/helpers.h"
# include "../arrays/matrix2D.h"
# include "../utils/concepts.h"

template <typename T1, typename T2>
void computeLU(std::vector< std::vector<T1> > & L, std::vector< std::vector<T1> > & U, 
        std::unique_ptr< std::vector< std::vector<T2> > > & myMatrix, std::size_t size) {
    std::size_t i{}, j{}, k{};
    typedef typename std::common_type<T1, T2>::type T;
    T sum{T(0)};
    for (i = 0; i < size; i++) {
        for (k = i; k < size; k++) {
            sum = T(0);
            for (j = 0; j < i; j++)
                sum += L.at(i).at(j)*U.at(j).at(k);
            U.at(i).at(k) = T(myMatrix -> at(i).at(k)) - sum;
         }
        for (k = i; k < size; k++) {
            if (i == k)
                L.at(i).at(i) = T(1);
            else {
                sum = T(0);
                for (j = 0; j < i; j++)
                    sum += L.at(k).at(j)*U.at(j).at(i);
                L.at(k).at(i) = (T(myMatrix -> at(k).at(i)) - sum)/U.at(i).at(i);
            }
        }        
    }   
}

template <typename T>
auto doolittleLUDecomposition(Matrix2D<T> & matrix) {
    std::size_t size = std::get<0>(matrix.getShape());
    if constexpr (std::is_integral<T>::value) {
        double sum{};
        std::unique_ptr<std::vector< std::vector<T>>> myMatrix = matrix.getCopyPtrMatrix();
        std::vector< std::vector<double>> L, U;
        U.resize(size, std::vector<double>(size, 0.0));
        L.resize(size, std::vector<double>(size, 0.0));
        computeLU(L, U, myMatrix, size);
        Matrix2D<double> upperTriang(U);
        Matrix2D<double> lowerTriang(L);
        return std::make_tuple(upperTriang, lowerTriang);        
    };
    if constexpr(std::is_floating_point<T>::value) {
        T sum{};
        std::unique_ptr<std::vector< std::vector<T>>> myMatrix = matrix.getCopyPtrMatrix();
        std::vector< std::vector<T>> L, U;
        U.resize(size, std::vector<T>(size, T(0)));
        L.resize(size, std::vector<T>(size, T(0)));
        computeLU(L, U, myMatrix, size);        
        Matrix2D<T> upperTriang(U);
        Matrix2D<T> lowerTriang(L);
        return std::make_tuple(upperTriang, lowerTriang);
    };
}

// Gauss Elimination for Ax = b solution //

template <typename T>
requires std::is_floating_point<T>::value && TwoD<Matrix2D<T>> && OneD<Matrix1D<T>>
std::pair<Matrix2D<T>, Matrix1D<T>> gaussianElimination(Matrix2D<T> A, Matrix1D<T> b) {
    std::size_t size = std::get<0>(A.getShape());
    std::size_t pivotIndex;
    double iRC{0.0};
    for (std::size_t k = 0; k < size; ++k) {
        pivotIndex = A.maxInColIndex(k);
        A.swapRows(pivotIndex, k);
        b.swapRows(pivotIndex, k);
        for (std::size_t i = k + 1; i < size; ++i) {
            iRC = double(A(i,k))/double(A(k,k));
            for (std::size_t j = k; j < size; ++j)
                A(i,j) = A(i,j) - iRC*A(k,j);
            b(i) = b(i) - iRC*b(i);
        }
    } 
    b(size - 1) = b(size - 1)/A(size - 1,size - 1);
    for (std::size_t k = size - 2; k > 0; k--) {
        for (std::size_t j = k + 1; j < size; j++)
            b(k) -= A(k, j)*b(j);
        b(k) = b(k)/A(k,k);
    }

    return std::make_pair(A, b);
}

# endif 