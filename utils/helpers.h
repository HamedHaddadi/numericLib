# ifndef HELPERS_H_
# define HELPERS_H_
# include "./header.h"
# include "./concepts.h"

// array multiplier
// works for any container that supports size()
// helper works with containers of same types

template <typename CONT>
void matrixMultiply(const CONT & matOne, const CONT & matTwo, CONT & results) {

    constexpr bool has2DSize = requires(const CONT & matrix) {
        matrix.size();
        matrix.at(0).size();
    };
    if constexpr (!has2DSize) 
        throw std::invalid_argument("the passed container does not support multiplication!");
    else {
        std::size_t row = matOne.size();
        std::size_t col = matTwo.at(0).size();
        std::size_t comOne = matOne.at(0).size();
        std::size_t comTwo = matTwo.size();
        if (comOne != comTwo)
            throw std::length_error("matrix length mismatch in matrixMultiply");

        for (std::size_t i = 0; i < row; ++i)
            for (std::size_t j = 0; j < col; ++j)
                for (std::size_t k = 0; k < comOne; ++k)
                    results.at(i).at(j) += matOne.at(i).at(k)*matTwo.at(k).at(j);
    }
}


template <typename T>
std::vector< std::vector<T>> operator+(const std::vector< std::vector<T>> & matOne, const std::vector< std::vector<T>> & matTwo) {

    std::size_t rowOne = matOne.size();
    std::size_t colOne = matOne.at(0).size();  
    std::size_t rowTwo = matTwo.size();
    std::size_t colTwo = matTwo.at(0).size();
    std::vector< std::vector<T>> results;
    results.resize(rowOne, std::vector<T>(colOne, T()));

    if ((rowOne != rowTwo) && (colOne != colTwo))
        throw std::length_error("matrix length mismatch in addition");
    else {
        for (std::size_t i = 0; i < rowOne; ++i)
            for (std::size_t j = 0; j < colOne; ++j)
                results.at(i).at(j) = matOne.at(i).at(j) + matTwo.at(i).at(j);
    }
    return results;
} 

template <typename T>
std::vector< std::vector<T> > operator-(const std::vector< std::vector<T>> & matOne, const std::vector< std::vector<T>> & matTwo) {

    std::size_t rowOne = matOne.size();
    std::size_t colOne = matOne.at(0).size();  
    std::size_t rowTwo = matTwo.size();
    std::size_t colTwo = matTwo.at(0).size();
    std::vector< std::vector<T>> results;
    results.resize(rowOne, std::vector<T>(colOne, T()));

    if ((rowOne != rowTwo) && (colOne != colTwo))
        throw std::length_error("matrix length mismatch in addition");
    else {
        for (std::size_t i = 0; i < rowOne; ++i)
            for (std::size_t j = 0; j < colOne; ++j)
                results.at(i).at(j) = matOne.at(i).at(j) - matTwo.at(i).at(j);
    }
    return results;
} 


template <typename T>
void lowerInit(std::vector< std::vector<T> > & matrix, std::size_t row, std::size_t col) {
    for (std::size_t iDx = 0; iDx < row; ++iDx)
        for (std::size_t iDy = 0; iDy < col; ++iDy) {
            if (iDx == iDy)
                matrix.at(iDx).at(iDy) = T(1);
            else if (iDy > iDx)
                matrix.at(iDx).at(iDy) = T(0);
        }
}

template <int I, typename ...Types>
decltype(auto) getFromTuple(Types&&... elements) {
    return std::get<I>(std::forward_as_tuple(elements...));
}

// useful methods for binary comparisons
template <typename T>
bool binaryComparison(T elemOne, T elemTwo) {

    if constexpr (std::is_integral<T>::value) {
        if (elemOne == elemTwo)
            return true;
        else
            return false;
    }
    else if constexpr (std::is_floating_point<T>::value) {
        T epsilon_{T(1e-6)};
        if (std::fabs(elemOne - elemTwo) <= epsilon_)
            return true;
        else
            return false;
    }

}

template <typename U, typename W>
bool binaryComparison(U elemOne, W elemTwo) {
    double epsilon_{1e-6};
    if (std::fabs(static_cast<double>(elemOne) - static_cast<double>(elemTwo)) <= epsilon_)
        return true;
    else    
        return false;
}

// check to see the ccontainers are square vectors
template <typename T>
requires isVector<T>
bool isSquare(T vec) {
    if (vec.size() == vec.at(0).size())
        return true;
    else
        return false;
}


class LibException: public std::exception {
    private:
        const char * message_;
    public:
        LibException(const char * msg):message_{msg}{}
        const char* what() const noexcept override {
            return message_;
        }
};








# endif 