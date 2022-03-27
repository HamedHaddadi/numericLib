
# ifndef ARRAYND_H_
# define ARRAYND_H_
# include "../utils/header.h"

template <typename T, int DIM> class Array;

/* + */
template <typename T, int DIM>
Array<T, DIM> operator+(const Array<T, DIM> &, const Array<T, DIM> &);

// + for chaining
template <typename T, int DIM>
Array<T, DIM>&& operator+(const Array<T, DIM> &, Array<T, DIM> &&);

template <typename T, int DIM>
Array<T, DIM>&& operator+(Array<T, DIM> &&, Array<T, DIM> &&);

/* - */
template <typename T, int DIM>
Array<T, DIM> operator-(const Array<T, DIM> &, const Array<T, DIM> &);

template <typename T, int DIM>
Array<T, DIM>&& operator-(const Array<T, DIM> &, Array<T, DIM> &&);

template <typename T, int DIM>
Array<T, DIM>&& operator-(Array<T, DIM> &&, Array<T, DIM> &&);

// *
template <typename U, typename V, int D>
Array<std::common_type_t<U,V>, D> operator*(const Array<U, D> &, const Array<V, D> &);

// /
template <typename T, int DIM>
auto operator/(const Array<T, DIM> &, const T);

// / overload
template <typename T, int DIM>
auto operator/(const Array<T, DIM> &, const Array<T, DIM> &);

/* << */
template <typename T, int DIM>
std::ostream& operator<<(std::ostream & , const Array<T, DIM> &); 

template <typename T, int DIM>
class Array {
    public:
        // type aliases
        using vector1d = std::vector<T>;
        using vec_it = typename std::vector<T>::iterator;
        using vec_cit = typename std::vector<T>::const_iterator;
        Array() = default;
        Array(auto row):row_{static_cast<std::size_t>(row)}{
            array_.resize(row_, T());
            size_ = row_;
        };
        Array(auto row, auto column):row_{static_cast<std::size_t>(row)}, column_{static_cast<std::size_t>(column)}{
            array_.resize(row*column, T());
            size_ = row_*column_;
        };
        Array(auto row, auto column, vector1d array):row_{static_cast<std::size_t>(row)},
             column_{static_cast<std::size_t>(column)}, array_{std::move(array)}{size_ = row_*column_;};
        // overloaded operators
        T& operator()(const auto, const auto);
        const T & operator()(const auto, const auto) const; 
        friend Array<T, DIM> operator+ <T> (const Array<T, DIM> &, const Array<T, DIM> &);
        friend Array<T, DIM>&& operator+ <T> (const Array<T, DIM> &, Array<T, DIM> &&);
        friend Array<T, DIM>&& operator+ <T> (Array<T, DIM> &&, Array<T, DIM> &&);        
        friend Array<T, DIM> operator- <T> (const Array<T, DIM> &, const Array<T, DIM> &);
        friend Array<T, DIM>&& operator- <T> (const Array<T, DIM> &, Array<T, DIM> &&);
        friend Array<T, DIM>&& operator- <T> (Array<T, DIM> &&, Array<T, DIM> &&);
        template <typename U, typename V, int D>
        friend Array<std::common_type_t<U,V>, D> operator* (const Array<U, D> &, const Array<V, D> &);
        friend auto operator/ <T> (const Array<T, DIM> &, const T);
        friend auto operator/ <T> (const Array<T, DIM> &, const Array<T, DIM> &);
        friend std::ostream & operator<< <T> (std::ostream &, const Array<T, DIM> &);
        // data access
        std::tuple<std::size_t, std::size_t> getSize() {return std::make_tuple(row_, column_);}
        // utilities
        void fillRandom(T, T);
        void inplaceTranspose();  
        Array<T, DIM> transpose(); 
        T det();
        Array<T, DIM> adjoint();
        auto inverse();
        static T computeDeterminant(const std::vector<T> &, std::size_t);
        static std::unique_ptr<std::vector<T>> coFactor(const std::vector<T> &, std::size_t, std::size_t, std::size_t);
        std::unique_ptr<std::vector< std::vector<T>>> getMatrixForm();
        // conditions
        bool isSingular();
        bool isSparse();

    private:
        int dim_{DIM};
        std::size_t row_{0}, column_{0};
        std::size_t size_{0};
        std::vector<T> array_;    
};

//operators
template <typename T, int DIM>
T& Array<T, DIM>::operator()(const auto x, const auto y) {
    return array_.at(static_cast<std::size_t>(x)*column_ + static_cast<std::size_t>(y));
}

template <typename T, int DIM>
const T& Array<T, DIM>::operator()(const auto x, const auto y) const {
    return array_.at(static_cast<std::size_t>(x)*column_ + static_cast<std::size_t>(y));
}

template <typename T, int DIM>
Array<T, DIM> operator+ (const Array<T, DIM> & arrOne, const Array<T, DIM> & arrTwo) {

    if (arrOne.size_ != arrTwo.size_) {
        throw std::length_error("array length mismatch for +");
    }
    Array<T, DIM> results(arrOne.row_, arrOne.column_);
    for (std::size_t s = 0; s < arrOne.size_; ++s)
        results.array_.at(s) = arrOne.array_.at(s) + arrTwo.array_.at(s);
    return results;
}

template <typename T, int DIM>
Array<T, DIM>&& operator+ (const Array<T, DIM> & arrOne, Array<T, DIM> && arrTwo) {

    if (arrOne.size_ != arrTwo.size_) {
        throw std::length_error("array length mismatch for +");
    }
    for (std::size_t s = 0; s < arrOne.size_; ++s)
        arrTwo.array_.at(s) += arrOne.array_.at(s);
    return std::move(arrTwo);
}

template <typename T, int DIM>
Array<T, DIM>&& operator+ (Array<T, DIM> && arrOne, Array<T, DIM> && arrTwo) {

    if (arrOne.size_ != arrTwo.size_) {
        throw std::length_error("array length mismatch for +");
    }
    for (std::size_t s = 0; s < arrOne.size_; ++s)
        arrOne.array_.at(s) += arrTwo.array_.at(s);
    return std::move(arrOne);
}

template <typename T, int DIM>
Array<T, DIM> operator - (const Array<T, DIM> & arrOne, const Array<T, DIM> & arrTwo) {
    if (arrOne.size_ != arrTwo.size_) {
        throw std::length_error("array length mismatch for +");
    }
    Array<T, DIM> results(arrOne.row_, arrOne.column_);
    for (std::size_t s = 0; s < arrOne.size_; ++s)
        results.array_.at(s) = arrOne.array_.at(s) - arrTwo.array_.at(s);
    return results;
}

template <typename T, int DIM>
Array<T, DIM>&& operator - (const Array<T, DIM> & arrOne, Array<T, DIM> && arrTwo) {
    if (arrOne.size_ != arrTwo.size_) {
        throw std::length_error("array length mismatch for +");
    }
    for (std::size_t s = 0; s < arrTwo.size_; ++s)
        arrTwo.array_.at(s) -= arrOne.array_.at(s);
    return std::move(arrTwo);
}

template <typename T, int DIM>
Array<T, DIM>&& operator - (Array<T, DIM> && arrOne, Array<T, DIM> && arrTwo) {
    if (arrOne.size_ != arrTwo.size_) {
        throw std::length_error("array length mismatch for +");
    }
    for (std::size_t s = 0; s < arrTwo.size_; ++s)
        arrTwo.array_.at(s) -= arrOne.array_.at(s);
    return std::move(arrTwo);
}

template <typename U, typename V, int D>
Array<std::common_type_t<U,V>, D> operator*(const Array<U, D> & arrOne, const Array<V, D> & arrTwo) {

    if (arrOne.column_ != arrTwo.row_) {
        throw std::length_error("array length mismatch for *");
    }

    Array<std::common_type_t<U,V>, D> results(arrOne.row_, arrTwo.column_);
    std::size_t indexResults{0}; 
    typedef typename std::common_type<U,V>::type T;
    T sum{T()};
    for (std::size_t x = 0; x < arrOne.row_; ++x) 
        for (std::size_t y = 0; y < arrTwo.column_; ++y) {
            sum  = T(0.0);
            for (std::size_t c = 0; c < arrOne.column_; ++c) {
                sum += arrOne.array_.at(arrOne.column_*x + c)*arrTwo.array_.at(arrTwo.column_*c + y);
        }
        results.array_.at(indexResults++) = sum;
    }
    return results;
}

template <typename T, int DIM>
auto operator/(const Array<T, DIM> & arr, const T denum) {
    std::vector<double> results(arr.size_, 0.0);
    float denumF = static_cast<double>(denum);
    for (std::size_t index = 0; index < arr.size_; ++index) 
        results.at(index) = static_cast<double>(arr.array_.at(index))/denum;
    Array<double, 2> resultArray(arr.row_, arr.column_, results);
    return resultArray;
}

template <typename T, int DIM>
auto operator/(const Array<T, DIM> & arrOne, const Array<T, DIM> & arrTwo) {
    std::vector<double> results(arrOne.size_);
    for (std::size_t index = 0; index <arrOne.size_; ++index)
        results.at(index) = static_cast<double>(arrOne.array_.at(index))/static_cast<double>(arrTwo.array_.at(index));
    Array<double, 2> resultArray(arrOne.row_, arrOne.columns_, results);    
    return resultArray;
}

template <typename T, int DIM>
std::ostream & operator<< (std::ostream & outStream, const Array<T, DIM> & arr) {
    for (auto s = 0; s < arr.size_; s++) {
        outStream<<arr.array_.at(s)<<" ";
        if ((arr.dim_ == 2) && ((s + 1) % (arr.column_) == 0) && (s > 0)) {
            outStream<<std::endl;
        }
        if ((arr.dim_ == 3) && (s + 1 % (arr.row_)*(arr.column_) == 0) && (s > 0))
            outStream<<"--------------"<<std::endl;
    }
    return outStream;
}

template <typename T, int DIM>
void Array<T, DIM>::fillRandom(T min, T max) {
    std::random_device seeder;
    std::mt19937 engine(seeder());
    if constexpr (std::is_integral<T>::value) {
        std::uniform_int_distribution<T> dist(min, max);
        for (auto s = 0; s < size_; ++s)
            array_.at(s) = dist(engine);
    }
    else if constexpr (std::is_floating_point<T>::value) {
        std::uniform_real_distribution<T> dist(min, max);
        for (auto s = 0; s < size_; ++s)
            array_.at(s) = dist(engine);
    }         
}

template <typename T, int DIM>
void Array<T, DIM>::inplaceTranspose() {
    std::vector<T> tempArray(size_, T());
    std::size_t index{0};
    for (std::size_t col = 0; col < column_; ++col)
        for (std::size_t row = 0; row < row_; ++row)
            tempArray.at(index++) = array_.at(row*column_ + col);
    std::swap(row_, column_);
    std::swap(array_, tempArray);
}

template <typename T, int DIM>
Array<T, DIM> Array<T, DIM>::transpose() {
    Array<T, DIM> transArray = *this;
    transArray.inplaceTranspose();
    return transArray;
}

// determinant calculation methods
template <typename T, int DIM>
std::unique_ptr<std::vector<T>> Array<T, DIM>::coFactor(const std::vector<T> & myArray, std::size_t currentRow, 
            std::size_t currentCol, std::size_t currentSize) {
    
    std::unique_ptr<std::vector<T>> coFactArray = std::make_unique<std::vector<T>>();
    coFactArray -> resize(myArray.size());
    std::size_t index{0};
    for (std::size_t idx = 0; idx < currentSize; ++idx)
        for (std::size_t idy = 0; idy < currentSize; ++idy) {
            if (idx != currentRow && idy != currentCol) {
                coFactArray -> at(index++) = myArray.at(idx*currentSize + idy);
            }
        }
    return std::move(coFactArray);
}

template <typename T, int DIM>
T Array<T, DIM>::computeDeterminant(const std::vector<T> & myArray, std::size_t colSize) {
    T detVal{};
    if (colSize == 1)
        return myArray.at(0);
    T sign{1};
    for (std::size_t col = 0; col < colSize; ++col) {
        detVal += sign*computeDeterminant(*coFactor(myArray, 0, col, colSize), colSize - 1)*myArray.at(col);
        sign = -sign;
    }
    return detVal;
}

template <typename T, int DIM>
T Array<T, DIM>::det() {
    TIME_THIS();
    T detValue{};
    detValue = computeDeterminant(array_, column_);
    return detValue;
}

template <typename T, int DIM>
Array<T, DIM> Array<T, DIM>::adjoint() {
    if (row_ != column_) {
        throw std::length_error("array must be square for finding adjoint matrix");
    }
    std::vector<T> tempVec(size_, T());
    int sign{1};
    for (std::size_t row = 0; row < row_; ++row)
        for (std::size_t col = 0; col < column_; ++col) {
            sign = ((row + col)%2 == 0) ? 1:-1;
            tempVec.at(row*column_ + col) = sign*computeDeterminant(*coFactor(array_, row, col, column_), column_ - 1);
        }
    Array<T, DIM> adjArray(row_, column_, tempVec);
    adjArray.inplaceTranspose();
    return adjArray;
}

template <typename T, int DIM>
auto Array<T, DIM>::inverse() {
    if (row_ != column_) {
        throw std::length_error("array must be square for finding inverse matrix");
    }
    T detValue = det();
    if (detValue == T(0)) {
        throw std::overflow_error("division by zero in determinant calculations");
    }

    Array<double, 2> invArray = adjoint()/detValue;
    return invArray;

}

template <typename T, int DIM>
bool Array<T, DIM>::isSingular() {
    T detValue{0};
    detValue = det();
    if constexpr (std::is_integral<T>::value) {
        if (detValue == 0)
            return true;
        else
            return false;
    }
    else {
        if (detValue < 2.0*std::numeric_limits<T>::min())
            return true;
        else
            return false;
    }
}

template <typename T, int DIM>
bool Array<T, DIM>::isSparse() {
    auto nonZero =  std::count_if(array_.begin(), array_.end(),[](T elem){ return elem > T(0);});
    if (static_cast<float>(nonZero) > static_cast<float>(0.5*size_))
        return false;
    else 
        return true;
}


template <typename T, int DIM>
std::unique_ptr<std::vector<std::vector<T>>> Array<T, DIM>::getMatrixForm() {
    using vec2D = std::vector< std::vector<T> >;
    std::unique_ptr<vec2D> matrix = std::make_unique<vec2D>();
    
    if constexpr (DIM == 2) {
        matrix -> resize(row_, std::vector<T>(column_, T()));
        std::size_t idx{0}, idy{0};
        for (std::size_t index = 0; index < size_; ++index) {
            matrix -> at(idx).at(idy++) = array_.at(index);
            if ((index + 1) % column_ == 0) {
                idy = 0;
                ++idx;
            }
        }
        return std::move(matrix);
    }
    else if constexpr (DIM == 1) {
        matrix -> resize(1, std::vector<T>(size_, T()));
        for (std::size_t index = 0; index < size_; ++index) {
            matrix -> at(0).at(index) = array_.at(index);
        }
        return std::move(matrix);
    }
}


# endif 