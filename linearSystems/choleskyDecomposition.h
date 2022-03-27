# ifndef CHOLESKYDECOMPOSITION_H_
# define CHOLESKYDECOMPOSITION_H_

# include "../utils/header.h"

//
template <typename T>
T getLSum(Matrix2D<T> & L, const std::size_t & iDy) {
    T tmpSum{0};
    for (std::size_t iD = 0; iD < iDy; ++iD)
        tmpSum += ::pow(L(iDy, iD), 2);
    return tmpSum;
}

template <typename T>
T getLLSum(Matrix2D<T> & L, const std::size_t & iDx, const std::size_t & iDy) {
    T tmpSum{0};
    for (std::size_t iD = 0; iD < iDy; ++iD)
        tmpSum += L(iDx, iD)*L(iDy, iD);
    return tmpSum;
}

template <typename T>
bool computeCholesky(Matrix2D<T> & L, Matrix2D<T> & A, const std::size_t & size) {
    T sum{0.0};
    for (std::size_t iDx = 0; iDx < size; ++iDx) {
        for (std::size_t iDy = 0; iDy <= iDx; ++iDy) {
            if (iDx == iDy) {
                sum = getLSum(L, iDy);
                if (A(iDy, iDy) >= sum)
                    L(iDy, iDy) = ::sqrt(A(iDy, iDy) - sum);
                else
                    return false;
            }
            else {
                sum = getLLSum(L, iDx, iDy);
                L(iDx, iDy) = (A(iDx, iDy) - sum)/L(iDy, iDy);
            }
        }
    }
    return true;
}

template <typename T>
auto choleskyDecomposition(Matrix2D<T> & A) {
    TIME_THIS();
    std::size_t size = std::get<0>(A.getShape());
    bool isPositiveDefinite{true};
    if constexpr (std::is_integral<T>::value) {
        Matrix2D<double> L(size, size);
        isPositiveDefinite = computeCholesky(L, A, size);
        return std::make_pair(isPositiveDefinite, L);
    }
    else if constexpr (std::is_floating_point<T>::value) {
        Matrix2D<T> L(size, size);
        isPositiveDefinite = computeCholesky(L, A, size);
        return std::make_pair(isPositiveDefinite, L);
    }

}

# endif 
