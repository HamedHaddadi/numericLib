# ifndef CONCEPTS_H_
# define CONCEPTS_H_
# include "../arrays/matrix2D.h"


template <typename U, typename W>
concept areSameDataType = std::is_same<typename U::value_type, typename W::value_type>::value;

template <typename U, typename W>
concept areSameVectors = (std::is_same<U, std::vector<int>>::value && 
                            std::is_same<W, std::vector<int>>::value) || 
                            (std::is_same<U, std::vector<double>>::value &&
                                std::is_same<W, std::vector<double>>::value);

template <typename T>
concept isVector = std::is_same<T, std::vector<int>>::value || std::is_same<T, std::vector<double>>::value;


template <typename Matrix>
concept TwoD = requires (Matrix mat) {
    Matrix::Dim == 2;
};

template <typename Matrix>
concept OneD = requires (Matrix mat) {
    Matrix::Dim == 1;
};


# endif 