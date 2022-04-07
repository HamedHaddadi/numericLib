# ifndef MATRIX1D_HH_
# define MATRIX1D_HH_

# include "../utils/header.h"

template <typename T> class Matrix1D;

template <typename T>
std::ostream& operator<<(std::ostream &, const Matrix1D<T> &);

template <typename T>
T dot(const Matrix1D<T> &, const Matrix1D<T> &);

template <typename T>
Matrix1D<T> operator+(const Matrix1D<T> &, const Matrix1D<T> &);

template <typename T>
Matrix1D<T> operator-(const Matrix1D<T> &, const Matrix1D<T> &);


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
        void fillRandom(T, T);
        void swapRows(std::size_t, std::size_t);
        std::size_t getShape() {return size_;}
        friend std::ostream & operator<< <T> (std::ostream &, const Matrix1D<T> &);
        // math operator
        friend T dot <T> (const Matrix1D<T> &, const Matrix1D<T> &);
        friend Matrix1D<T> operator+ <T> (const Matrix1D<T> &, const Matrix1D<T> &);
        friend Matrix1D<T> operator- <T> (const Matrix1D<T> &, const Matrix1D<T> &);
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
T dot(const Matrix2D<T> & matOne, const Matrix1D<T> & matTwo) {
    T init{T(0)};
    std::inner_product(matOne.matrix_.begin(), matOne.matrix_.end(), matTwo.matrix_.begin(), init);
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


# endif 