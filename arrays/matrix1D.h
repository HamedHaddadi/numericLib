# ifndef MATRIX1D_HH_
# define MATRIX1D_HH_

# include "../utils/header.h"
# include "../utils/concepts.h"

template <typename T> class Matrix1D;

template <typename T>
std::ostream& operator<<(std::ostream &, const Matrix1D<T> &);

template <typename T>
T dot(const Matrix1D<T> &, const Matrix1D<T> &);

template <typename T>
Matrix1D<T> operator+(const Matrix1D<T> &, const Matrix1D<T> &);

template <typename T>
Matrix1D<T> operator-(const Matrix1D<T> &, const Matrix1D<T> &);

// norms
template <typename T>
T normInfty(const Matrix1D<T> &);

template <typename T>
T norm1(const Matrix1D<T> &);

template <typename T>
T norm2(const Matrix1D<T> &);

template <typename T>
T euclidean(const Matrix1D<T> &, const Matrix1D<T> &);

template <typename T>
bool sameShape(const Matrix1D<T> &, const Matrix1D<T> &);

template <typename U, typename W>
bool sameShape(const Matrix1D<U> &, const Matrix1D<W> &);

template <typename U, typename W>
bool operator==(const Matrix1D<U> &, const Matrix1D<W> &);

template <typename T>
class Matrix1D {
    private:
        std::size_t size_{0};
        std::vector<T> matrix_;
    public:
        using matrixType = T;
        inline static constexpr int Dim{1};
        inline static constexpr bool FixedSize{false}; 
        inline static constexpr bool linearized{true};        
        Matrix1D() = default;
        Matrix1D(std::size_t);
        Matrix1D(std::vector<T>);
        ~Matrix1D() = default;
        T& operator() (std::size_t);
        const T& operator()(std::size_t) const;
        void fillRandom(T, T);
        void swapRows(std::size_t, std::size_t);
        void inputFromString(std::string);
        T norm1();
        T norm2();
        T normInfty();
        std::size_t getShape() {return size_;}
        const std::size_t getShape() const {return size_;}
        friend std::ostream & operator<< <T> (std::ostream &, const Matrix1D<T> &);
        // math operator
        friend T dot <T> (const Matrix1D<T> &, const Matrix1D<T> &);
        friend Matrix1D<T> operator+ <T> (const Matrix1D<T> &, const Matrix1D<T> &);
        friend Matrix1D<T> operator- <T> (const Matrix1D<T> &, const Matrix1D<T> &);
        friend bool sameShape <T> (const Matrix1D<T> &, const Matrix1D<T> &);
        template <typename U, typename W>
        friend bool sameShape (const Matrix1D<U> &, const Matrix1D<W> &);
        template <typename U, typename W>
        friend bool operator==(const Matrix1D<U> &, const Matrix1D<W> &);
        // Norms
        friend T normInfty <T> (const Matrix1D<T> &);
        friend T norm1 <T> (const Matrix1D<T> &);
        friend T norm2 <T> (const Matrix1D<T> &);
        friend T euclidean <T> (const Matrix1D<T> &, const Matrix1D<T> &);
};

template <typename T>
Matrix1D<T>::Matrix1D(std::size_t size):size_{size} {
    matrix_.resize(size_, T(0));
}

template <typename T>
Matrix1D<T>::Matrix1D(std::vector<T> vec):matrix_{std::move(vec)} {
    size_ = vec.size();
}

template <typename T>
T& Matrix1D<T>::operator() (std::size_t index) {
    return matrix_.at(index);
}

template <typename T>
const T& Matrix1D<T>::operator() (std::size_t index) const {
    return matrix_.at(index);
}

template <typename T>
void Matrix1D<T>::fillRandom(T min, T max) {
    std::random_device seeder;
    std::mt19937 engine(seeder());
    if constexpr (std::is_integral<T>::value) {
        std::uniform_int_distribution<T> dist(min, max);
        for (std::size_t i = 0; i < size_; i++)
            matrix_.at(i) = dist(engine);
    }
    else if constexpr (std::is_floating_point<T>::value) {
        std::uniform_real_distribution<T> dist(min, max);
        for (std::size_t i = 0; i < size_; i++)
            matrix_.at(i) = dist(engine);
    }               
}

// norms
template <typename T>
T Matrix1D<T>::norm1() {
    T normValue{T()};
    std::for_each(matrix_.begin(), matrix_.end(), [&](T elem) {
        normValue += std::abs(elem);
    });
    return normValue;
}

template <typename T>
T Matrix1D<T>::norm2() {
    T normValue{T()};
    normValue = std::inner_product(matrix_.begin(), matrix_.end(), matrix_.begin(), normValue);
    return ::sqrt(normValue);
}

template <typename T>
T Matrix1D<T>::normInfty() {
    T maxValue = std::numeric_limits<T>::min();
    std::for_each(matrix_.begin(), matrix_.end(), [&](T elem){
        if (std::abs(elem) >= maxValue) {
            maxValue = std::abs(elem);
        }
    });
    return maxValue;
}

template <typename T>
void Matrix1D<T>::swapRows(std::size_t rowOne, std::size_t rowTwo) {
    std::swap(matrix_.at(rowOne), matrix_.at(rowTwo));
}

template <typename T>
std::ostream& operator<<(std::ostream & outStream, const Matrix1D<T> & mat) {

    for (std::size_t size = 0; size < mat.size_; ++size)
        outStream<< mat.matrix_.at(size)<<" ";
    outStream<<std::endl;
    return outStream;

} 

template <typename T>
T dot(const Matrix1D<T> & matOne, const Matrix1D<T> & matTwo) {
    T init{T(0)};
    init = std::inner_product(matOne.matrix_.begin(), matOne.matrix_.end(), matTwo.matrix_.begin(), init);
    return init;
}

template <typename T>
Matrix1D<T> operator+(const Matrix1D<T> & matOne, const Matrix1D<T> & matTwo) {
    Matrix1D<T> sumResults(matOne.size_);
    std::transform(matOne.matrix_.begin(), matOne.matrix_.end(), matTwo.matrix_.begin(), sumResults.matrix_.begin(), std::plus<T>());
    return sumResults;
}

template <typename T>
Matrix1D<T> operator-(const Matrix1D<T> & matOne, const Matrix1D<T> & matTwo) {
    Matrix1D<T> minusResults(matOne.size_);
    std::transform(matOne.matrix_.begin(), matOne.matrix_.end(), matTwo.matrix_.begin(), minusResults.matrix_.begin(), std::minus<T>());
    return minusResults;    
}

template <typename T>
T normInfty(const Matrix1D<T> & matrix) {
    T maxValue = std::numeric_limits<T>::min();
    std::for_each(matrix.begin(), matrix.end(), [&](T elem) {
        if (std::abs(elem) >= maxValue) {
            maxValue = std::abs(elem);
        }
    });
    return maxValue;
}

template <typename T>
T norm1(const Matrix1D<T> & matrix) {
    T normValue{T()};
    normValue = std::for_each(matrix.begin(), matrix.end(), [&](T elem) {normValue += std::abs(elem);});
    return normValue;
}

template <typename T>
T norm2(const Matrix1D<T> & matOne) {
    T normValue{T()};
    normValue = dot(matOne, matOne);
    return ::sqrt(normValue);
}

template <typename T>
T euclidean(const Matrix1D<T> & matOne, const Matrix1D<T> & matTwo) {
    Matrix1D<T> distMat(matOne.size_);
    T dotProd{T(0)};
    distMat = matOne - matTwo;
    dotProd = dot(distMat, distMat);
    return ::sqrt(dotProd);
}

template <typename T>
void Matrix1D<T>::inputFromString(std::string inputString) {
    T inDigit{T()};
    std::size_t rowIdx{0};
    std::stringstream inStr(inputString);
    while (inStr >> inDigit) {
        matrix_.at(rowIdx++) = inDigit;
    }
}

template <typename T>
bool sameShape(const Matrix1D<T> & matOne, const Matrix1D<T> & matTwo) {
    if (matOne.size_ == matTwo.size_)
        return true;
    else
        return false;
}

template <typename U, typename W>
bool sameShape(const Matrix1D<U> & matOne, const Matrix1D<W> & matTwo) {
    if (matOne.size_ == matTwo.size_)
        return true;
    else
        return false;
}

template <typename U, typename W>
bool operator==(const Matrix1D<U> & matOne, const Matrix1D<W> & matTwo) {
    bool shapeEqual{true};
    bool elemEqual{false};
    if (!sameShape(matOne, matTwo)) {
        shapeEqual = false;
    }
    elemEqual = std::equal(matOne.matrix_.begin(), matOne.matrix_.end(), matTwo.matrix_.begin(), binaryComparison);
    if (elemEqual && shapeEqual)
        return true;
    else
        return false;
}

# endif 