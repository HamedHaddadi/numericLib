# ifndef MATRIX2D_H_
# define MATRIX2D_H_

# include "../utils/header.h"
# include "../utils/helpers.h"
// # include "../linearSystems/luSolvers.h"

// friend overloaded methods 
template<typename T> class Matrix2D;
// +
template<typename T>
Matrix2D<T> operator+(const Matrix2D<T> &, const Matrix2D<T> &);
template<typename T>
Matrix2D<T>&& operator+(const Matrix2D<T> &, Matrix2D<T> &&);
template<typename T>
Matrix2D<T>&& operator+(Matrix2D<T> &&, const Matrix2D<T> &);
template<typename T>
Matrix2D<T>&& operator+(Matrix2D<T> &&, Matrix2D<T> &&);
// -
template<typename T>
Matrix2D<T> operator-(const Matrix2D<T> &, const Matrix2D<T> &);
template<typename T>
Matrix2D<T>&& operator-(const Matrix2D<T> &, Matrix2D<T> &&);
template<typename T>
Matrix2D<T>&& operator-(Matrix2D<T> &&, const Matrix2D<T> &);
template<typename T>
Matrix2D<T>&& operator-(Matrix2D<T> &&, Matrix2D<T> &&);

// binary comparison operators
template <typename T>
bool operator==(const Matrix2D<T> &, const Matrix2D<T> &);

template <typename U, typename W>
bool operator==(const Matrix2D<U> &, const Matrix2D<W> &);

template <typename T>
bool sameShape(const Matrix2D<T> &, const Matrix2D<T> &);

template <typename U, typename W>
bool sameShape(const Matrix2D<U> &, const Matrix2D<W> &);

// x using Strassen for square (N x N) matrices with N % 2 = 0
template<typename U, typename W>
Matrix2D<std::common_type_t<U, W>> operator*(const Matrix2D<U> &, const Matrix2D<W> &);

// helper friends to carry-out multiplication
template <typename U, typename W>
void matmul(const Matrix2D<U> &, const Matrix2D<W> &, Matrix2D<std::common_type_t<U, W>> &);

template <typename U, typename W>
void strassen(const Matrix2D<U> &, const Matrix2D<W> &, Matrix2D<std::common_type_t<U, W>> &);

// inputs; outputs
template <typename T>
std::ofstream& operator<<(std::ofstream &, const Matrix2D<T> &);
template <typename T>
std::ostream&  operator<<(std::ostream &, const Matrix2D<T> &);

// matrix norms
template <typename T>
T norm2(const Matrix2D<T> &);

template <typename T>
T euclidean(const Matrix2D<T> &, const Matrix2D<T> &);

// useful helper methods
template <typename T>
T sumElements(const Matrix2D<T> &);

// Class declarations
template <typename T>
class Matrix2D {

    protected:
        std::size_t row_{0}, col_{0};
        std::vector< std::vector<T> > matrix_;
        typedef typename std::vector< std::vector<T> >::iterator rowIter_;
        typedef typename std::vector<T>::iterator colIter_;    

    public:
        using matrixType = T;
        inline static constexpr int Dim{2};
        inline static constexpr bool FixedSize{false}; 
        inline static constexpr bool linearized{false};
        Matrix2D() = default;
        Matrix2D(std::size_t, std::size_t);
        Matrix2D(std::vector< std::vector<T>>);
        Matrix2D& operator=(const Matrix2D & src) = default;
        Matrix2D& operator=(Matrix2D && src) noexcept;
        Matrix2D(const Matrix2D<T> &) = default;
        Matrix2D(Matrix2D<T> &&) noexcept;
        virtual ~Matrix2D() = default;   
        // binary relational operators 
        friend bool operator== <T> (const Matrix2D<T> &, const Matrix2D<T> &);
        
        template <typename U, typename W>
        friend bool operator==(const Matrix2D<U> &, const Matrix2D<W> &);

        friend bool sameShape <T> (const Matrix2D<T> &, const Matrix2D<T> &);

        template <typename U, typename W>
        friend bool sameShape (const Matrix2D<U> &, const Matrix2D<W> &);

        // +
        friend Matrix2D<T> operator+ <T> (const Matrix2D<T> &, const Matrix2D<T> &);
        friend Matrix2D<T>&& operator+ <T> (const Matrix2D<T> &, Matrix2D<T> &&);
        friend Matrix2D<T>&& operator+ <T> (Matrix2D<T> &&, const Matrix2D<T> &);
        friend Matrix2D<T>&& operator+ <T> (Matrix2D<T> &&, Matrix2D<T> &&);
        // -
        friend Matrix2D<T> operator- <T> (const Matrix2D<T> &, const Matrix2D<T> &);
        friend Matrix2D<T>&& operator- <T> (const Matrix2D<T> &, Matrix2D<T> &&);
        friend Matrix2D<T>&& operator- <T> (Matrix2D<T> &&, const Matrix2D<T> &);
        friend Matrix2D<T>&& operator- <T> (Matrix2D<T> &&, Matrix2D<T> &&);
        // *
        template <typename U, typename W>
        friend Matrix2D<std::common_type_t<U, W>> operator* (const Matrix2D<U> &, const Matrix2D<W> &);

        template <typename U, typename W>
        friend void matmul(const Matrix2D<U> &, const Matrix2D<W> &, Matrix2D<std::common_type_t<U, W>> &);

        template <typename U, typename W>
        friend void strassen(const Matrix2D<U> &, const Matrix2D<W> &, Matrix2D<std::common_type_t<U, W>> &);

        // input/output friends
        friend std::ofstream& operator<< <T> (std::ofstream &, const Matrix2D<T> &);
        friend std::ostream& operator<< <T> (std::ostream &, const Matrix2D<T> &);  

        // norms
        friend T norm2 <T> (const Matrix2D<T> &);
        friend T euclidean <T> (const Matrix2D<T> &, const Matrix2D<T> &);
        friend T sumElements <T> (const Matrix2D<T> &);

        T& operator() (std::size_t i, std::size_t j); 
        const T& operator()(std::size_t, std::size_t) const;
        std::vector<T> values() const;
        const std::vector< std::vector<T> > & getConstMatrix() const {return matrix_;}  
        std::vector< std::vector<T> > & getMutableMatrix() {return matrix_;}
        std::unique_ptr< std::vector< std::vector<T>>> getCopyPtrMatrix();
        void setMatrix(std::vector< std::vector<T> > &&);
        void setMatrix(std::vector< std::vector<T> > &);
        Matrix2D<T> transpose();
        void inplaceTranspose();
        //  auto adjoint();
        T det(std::string flag = "recursive");
        virtual T sum();
        std::tuple<std::size_t, std::size_t> getShape() const {return std::make_pair(row_, col_);}
        bool isSquare() const;
        bool isSymmetric() const;
        bool isDiagonalDominant() const;
        void fill(T);
        void fillRandom(T, T);
        void inputFromString(std::string);
        // non-static utilities
        void swapRows(int, int);
        std::size_t maxInColIndex(std::size_t);
        // static methods 
        static T recursiveDeterminant(const std::vector<std::vector<T> > &, std::size_t);
        static std::unique_ptr<std::vector< std::vector<T> > > coFactor(const std::vector< std::vector<T>> &, 
            std::size_t, std::size_t, std::size_t);
                
};

template <typename T>
Matrix2D<T>::Matrix2D(std::size_t row, std::size_t col):row_{static_cast<std::size_t>(row)},
    col_{static_cast<std::size_t>(col)} {
    matrix_.resize(row, std::vector<T>(col, T()));
}

template <typename T>
Matrix2D<T>::Matrix2D(std::vector<std::vector<T>> vec):matrix_{std::move(vec)}{
    row_ = matrix_.size();
    col_ = matrix_.at(0).size();
}

template <typename T>
Matrix2D<T>::Matrix2D(Matrix2D<T> && src) noexcept {
    std::swap(row_, src.row_);
    std::swap(col_, src.col_);
    std::swap(matrix_, src.matrix_);
}

template <typename T>
Matrix2D<T>& Matrix2D<T>::operator=(Matrix2D<T> && src) noexcept {
    if (this != &src) {
        std::swap(row_, src.row_);
        std::swap(col_, src.col_);
        std::swap(matrix_, src.matrix_);
        return *this;
    }
} 

template <typename T>
void Matrix2D<T>::fill(T value) {
    for (auto & row: matrix_)
        for (auto & col: row)
            col = value;
}

template <typename T>
T& Matrix2D<T>::operator()(std::size_t i, std::size_t j) {
    return matrix_.at(i).at(j);
}

template <typename T>
const T& Matrix2D<T>::operator()(std::size_t i, std::size_t j) const {
    return matrix_.at(i).at(j);
}

// overloaded operators
template <typename T>
bool sameShape(const Matrix2D<T> & matOne, const Matrix2D<T> & matTwo) {
    if ((matOne.row_ == matTwo.row) && (matOne.col_ == matTwo.col_))
        return true;
    else
        return false;
}

template <typename U, typename W>
bool sameShape(const Matrix2D<U> & matOne, const Matrix2D<W> & matTwo) {
    if ((matOne.row_ == matTwo.row) && (matOne.col_ == matTwo.col_))
        return true;
    else
        return false;
}

template <typename T>
bool operator==(const Matrix2D<T> & leftMat, const Matrix2D<T> & rightMat) {

    bool vecEqual{true};
    if (!sameShape(leftMat, rightMat))
        return false;
    else {
    for (std::size_t row = 0; row < leftMat.ow_; ++row) {
        if (!std::equal(leftMat.matrix_.at(row).begin(), leftMat.matrix_.at(row).end(), rightMat.matrix_.at(row).begin(), binaryComparison)) {
            vecEqual = false;
            break;
            }
        }
    }
    return vecEqual;
}

template <typename U, typename W>
bool operator==(const Matrix2D<U> & leftMat, const Matrix2D<W> & rightMat) {
    
    bool vecEqual{true};
    if (!sameShape(leftMat, rightMat))
        return false;
    else {
     for (std::size_t row = 0; row < leftMat.row_; ++row) {
        if (!std::equal(leftMat.matrix_.at(row).begin(), leftMat.matrix_.at(row).end(), rightMat.matrix_.at(row).begin(), binaryComparison)) {
            vecEqual = false;
            break;
            }
        }       
    }
    return vecEqual;
}


template <typename T>
Matrix2D<T> operator+(const Matrix2D<T> & T1, const Matrix2D<T> & T2) {
    
    std::size_t row, col;
    if (T1.shape() != T2.shape())
        throw std::length_error("array size mismatch in & + & operation");
    else {
        row = T1.row_;
        col = T1.col_;
    }

    Matrix2D<T> new_mat(row, col);
    for (std::size_t i = 0; i < row; i++)
        for (std::size_t j = 0; j < col; j++)
            new_mat.matrix_.at(i).at(j) = T1.matrix_.at(i).at(j) + T2.matrix_.at(i).at(j);
    return new_mat;
}

template<typename T>
Matrix2D<T>&& operator+(const Matrix2D<T> & T1, Matrix2D<T> && T2) {

    std::size_t row, col;
    if (T1.shape() != T2.shape())
        throw std::length_error("array size mismatch in & + && operation");
    else {
        row = T1.row_;
        col = T1.col_;
    }

    for (std::size_t i = 0; i < row; ++i)
        for (std::size_t j = 0; j < col; ++j)
            T2.matrix_.at(i).at(j) += T1.matrix_.at(i).at(j);

    return std::move(T2);

}

template<typename T>
Matrix2D<T>&& operator+(Matrix2D<T> && T1, const Matrix2D<T> & T2) {

    std::size_t row, col;
    assert (T1.shape() != T2.shape());

    row = T1.row_;
    col = T1.col_;

    for (std::size_t i = 0; i < row; ++i)
        for (std::size_t j = 0; j < col; ++j)
            T1.matrix_.at(i).at(j) += T2.matrix_.at(i).at(j);
    
    return std::move(T1);

}

template<typename T>
Matrix2D<T>&& operator+(Matrix2D<T> && T1, Matrix2D<T> && T2) {

    std::size_t row, col;
    
    assert (T1.shape() == T2.shape());
    
    row = T1.row_;
    col = T1.col_;
    
    for (std::size_t i = 0; i < row; ++i)
        for (std::size_t j = 0; j < col; ++j)
            T1.matrix_.at(i).at(j) += T2.matrix_.at(i).at(j);
    
    return std::move(T1);
}

// -
template <typename T>
Matrix2D<T> operator-(const Matrix2D<T> & T1, const Matrix2D<T> & T2) {
    
    std::size_t row, col;
    assert (T1.shape() == T2.shape());
    
    row = T1.row_;
    col = T1.col_;
    
    Matrix2D<T> new_mat(row, col);
    for (std::size_t i = 0; i < row; i++)
        for (std::size_t j = 0; j < col; j++)
            new_mat.matrix_.at(i).at(j) = T1.matrix_.at(i).at(j) - T2.matrix_.at(i).at(j);
    return new_mat;
}

template<typename T>
Matrix2D<T>&& operator-(const Matrix2D<T> & T1, Matrix2D<T> && T2) {

    std::size_t row, col;
    
    assert (T1.shape() == T2.shape());

    row = T1.row_;
    col = T1.col_;
    
    for (std::size_t i = 0; i < row; ++i)
        for (std::size_t j = 0; j < col; ++j)
            T2.matrix_.at(i).at(j) -= T1.matrix_.at(i).at(j);

    return std::move(T2);

}

template<typename T>
Matrix2D<T>&& operator-(Matrix2D<T> && T1, const Matrix2D<T> & T2) {

    std::size_t row, col;

    assert (T1.shape() == T2.shape());
    row = T1.row_;
    col = T1.col_;
    
    for (std::size_t i = 0; i < row; ++i)
        for (std::size_t j = 0; j < col; ++j)
            T1.matrix_.at(i).at(j) -= T2.matrix_.at(i).at(j);
    
    return std::move(T1);

}

template<typename T>
Matrix2D<T>&& operator-(Matrix2D<T> && T1, Matrix2D<T> && T2) {

    std::size_t row, col;
    assert (T1.shape() == T2.shape());
    
    row = T1.row_;
    col = T1.col_;
    for (std::size_t i = 0; i < row; ++i)
        for (std::size_t j = 0; j < col; ++j)
            T1.matrix_.at(i).at(j) -= T2.matrix_.at(i).at(j);
    
    return std::move(T1);
}

// definition of all * related methods
template<typename U, typename W>
// regular matmul
void matmul(const Matrix2D<U> & matOne, const Matrix2D<W> & matTwo, Matrix2D<std::common_type_t<U, W>> & results) {
    std::size_t row1, col1, row2, col2;
    for (std::size_t i = 0; i < matOne.row_; ++i)
        for (std::size_t k = 0; k < matTwo.col_; ++k)
            for (std::size_t j = 0; j < matOne.col_; ++j)
                results.matrix_.at(i).at(k) += matOne.matrix_.at(i).at(j)*matTwo.matrix_.at(j).at(k);
}

// strassen
template <typename U, typename W>
void strassen(const Matrix2D<U> & matOne, const Matrix2D<W> & matTwo, Matrix2D<std::common_type_t<U, W>> & results) {
    typedef typename std::common_type_t<U,W> R;
    using vector2D = std::vector< std::vector<R>>;
    using vector1D = std::vector<R>;
    std::size_t row{matOne.row_}, col{matOne.col_};
    vector2D A11, A12, A21, A22, B11, B12, B21, B22, M1, M2, M3, M4, M5, M6, M7, C11, C12, C21, C22, C;
    
    A11.resize(row/2, vector1D(col/2, R()));
    A12.resize(row/2, vector1D(col/2, R()));
    A21.resize(row/2, vector1D(col/2, R()));
    A22.resize(row/2, vector1D(col/2, R()));
    B11.resize(row/2, vector1D(col/2, R()));
    B12.resize(row/2, vector1D(col/2, R()));
    B21.resize(row/2, vector1D(col/2, R()));
    B22.resize(row/2, vector1D(col/2, R()));

    for  (std::size_t i = 0; i < row/2; ++i)
        for (std::size_t j = 0; j < col/2; j++) {
            A11.at(i).at(j) = matOne.matrix_.at(i).at(j);
            A12.at(i).at(j) = matOne.matrix_.at(i).at(col/2 + j);
            A21.at(i).at(j) = matOne.matrix_.at(row/2 + i).at(j);
            A22.at(i).at(j) = matOne.matrix_.at(row/2 + i).at(col/2 + j);
            B11.at(i).at(j) = matTwo.matrix_.at(i).at(j);
            B12.at(i).at(j) = matTwo.matrix_.at(i).at(col/2 + j);
            B21.at(i).at(j) = matTwo.matrix_.at(row/2 + i).at(j);
            B22.at(i).at(j) = matTwo.matrix_.at(row/2 + i).at(col/2 + j);            
        } 
    
    matrixMultiply((A11 + A22), (B11 + B22), M1);
    matrixMultiply((A21 + A22), B11,M2);
    matrixMultiply(A11, (B12 - B22),M3);
    matrixMultiply(A22, (B21 - B11),M4);
    matrixMultiply((A11 + A12), B22, M5);
    matrixMultiply((A21 - A11), (B11 + B12), M6);
    matrixMultiply((A12 - A22), (B21 + B22), M7);

    C11 = M1 + M4 - M5 + M7;
    C12 = M3 + M5;
    C21 = M2 + M4;
    C22 = M1 - M2 + M3 + M6;

    C.resize(row, vector1D(col, R()));
    for (std::size_t i = 0; i < row/2; ++i)
        for (std::size_t j = 0; j < col/2; ++j) {
            C.at(i).at(j) = C11.at(i).at(j);
            C.at(i).at(col/2 + j) = C12.at(i).at(j);
            C.at(row/2 + i).at(j) = C21.at(i).at(j);
            C.at(row/2 + i).at(col/2 + j) = C22.at(i).at(j);
        } 
    results.matrix_ = C;
}

template<typename U, typename W>
Matrix2D<std::common_type_t<U, W>> operator*(const Matrix2D<U> & matOne, const Matrix2D<W> & matTwo) {
    Matrix2D<std::common_type_t<U,W>> results(matOne.row_, matTwo.col_);
    if (!(matOne.row_ == matOne.col_ == matTwo.row_ == matTwo.col_) || !(matOne.row_ % 2 == 0) ||
             !(matTwo.col_ % 2 == 0)) {
                matmul(matOne, matTwo, results);
             }
    else {
        strassen(matOne, matTwo, results);
    }
    return results;
}

// input/output operations
template<typename T>
std::ofstream& operator<<(std::ofstream & outStream, const Matrix2D<T> & T1) {

    for (std::size_t row = 0; row < T1.row_; ++row)
        for (std::size_t col = 0; col < T1.col_; ++col) {
            outStream << T1.matrix_.at(row).at(col) <<" ";
            if (col == T1.col_ - 1)
                outStream<<std::endl;
        }   
    return outStream;
}

template<typename T>
std::ostream& operator<<(std::ostream & outStream, const Matrix2D<T> & T1) {

    for (std::size_t row = 0; row < T1.row_; ++row)
        for (std::size_t col = 0; col < T1.col_; ++col) {
            outStream << T1.matrix_.at(row).at(col) <<" ";
            if (col == T1.col_ - 1)
                outStream<<std::endl;
        }
    return outStream;
}

// norms

// useful helpers
template <typename T>
T sumElements(const Matrix2D<T> & matrix) {
    T sumValue{T()};
    std::for_each(matrix.matrix_.begin(), matrix.matrix_.end(), [& sumValue](typename Matrix2D<T>::rowIter_ row) {
            std::for_each(row -> begin(), row -> end(), [&](typename Matrix2D<T>::colIter_ col) {
            sumValue += *col;
        });
    }); 
    return sumValue; 
}

template <typename T>
T norm2(const Matrix2D<T> & matrix) {
    T sumElems{T()};
    sumElems = sumElements(matrix);
    return ::sqrt(sumElems);
}

template <typename T>
T euclidean (const Matrix2D<T> & matOne, const Matrix2D<T> & matTwo) {
    Matrix2D<T> distMatrix(matOne.row_, matOne.col_);
    distMatrix = matOne - matTwo;
    T sumElems{T()};
    sumElems = sumElements(distMatrix);
    return ::sqrt(sumElems); 
}

// other class methods
template <typename T>
std::vector<T> Matrix2D<T>::values() const {
    std::vector<T> values; 
    for (auto const & elem: matrix_)
        values.insert(values.end(), elem.begin(), elem.end());
    return values;
}

template<typename T>
T Matrix2D<T>::sum() {
    T total = 0;
    for (const auto & row: matrix_)
        for (const auto & col: row)
            total += col;
    return total;

}

template<typename T>
Matrix2D<T> Matrix2D<T>::transpose() {
    Matrix2D<T> transpose(col_, row_);
    for (std::size_t row = 0; row < row_; ++row)
        for (std::size_t col = 0; col < col_; ++col)
            transpose.matrix_.at(col).at(row) = matrix_.at(row).at(col);
    return transpose;
}


template<typename T>
void Matrix2D<T>::inplaceTranspose() {
    std::vector< std::vector<T> > temp_matrix;
    std::size_t temp_col = col_;
    std::size_t temp_row = row_;
    temp_matrix.resize(col_, std::vector<T>(row_));
    for (std::size_t row = 0; row < row_; row++)
        for (std::size_t col = 0; col < col_; col++)
            temp_matrix.at(col).at(row) = matrix_.at(row).at(col);
    matrix_ = temp_matrix;
    row_ = temp_col;
    col_ = temp_row;

}

template<typename T>
bool Matrix2D<T>::isSquare() const {
    if (row_ == col_) 
        return true;
    else
        return false;
}

template <typename T>
bool Matrix2D<T>::isSymmetric() const {
    bool is_sym = true;
    for (std::size_t i = 0; i < row_; i++)
        for (std::size_t j = 0; j < col_; j++) 
            if (matrix_.at(i).at(j) != matrix_.at(j).at(i)) 
                is_sym = false;
    return is_sym;
}

template <typename T>
bool Matrix2D<T>::isDiagonalDominant() const {
    bool isDominant{true};
    T sumCol{T()};
    std::size_t row, col;
    for (row = 0; row < row_; ++row) {
        for (col = 0; col < col_; ++col) {
            sumCol += std::fabs(matrix_.at(row).at(col));
        }
        if (!(std::fabs(matrix_.at(row).at(col)) >= sumCol)) {
            isDominant = false;
            break;
        }
    }
    return isDominant;
}

template<typename T>
void Matrix2D<T>::fillRandom(T min, T max) {
    std::random_device seeder;
    std::mt19937 engine(seeder());
    if constexpr (std::is_integral<T>::value) {
        std::uniform_int_distribution<T> dist(min, max);
        for (std::size_t i = 0; i < row_; i++)
            for (std::size_t j = 0; j < col_; j++) {
                matrix_.at(i).at(j) = dist(engine);
            }
    }
    
    else if constexpr (std::is_floating_point<T>::value) {
        std::uniform_real_distribution<T> dist(min, max);
        for (std::size_t i = 0; i < row_; i++)
            for (std::size_t j = 0; j < col_; j++) {
                // auto val = dist(engine);
                // std::cout <<"the value is "<< val <<std::endl;
                matrix_.at(i).at(j) = dist(engine);
            }
    }
}

template <typename T>
void Matrix2D<T>::inputFromString(std::string inputString) {
    T inDigit{T()};
    std::size_t rowIdx{0}, colIdx{0};
    std::stringstream inStr(inputString); 
    while (inStr >> inDigit) {
        matrix_.at(rowIdx).at(colIdx) = inDigit;
        if (colIdx == col_ - 1) {
            ++rowIdx;
            colIdx = 0;
        }
        else {
            ++colIdx;
        }        
    }
}

// setMatrix methods
template <typename T>
void Matrix2D<T>::setMatrix(std::vector< std::vector<T> > && new_matrix) {
    std::size_t row, col;
    row = new_matrix.size();
    col = new_matrix[0].size();
    row_ = row;
    col_ = col;
    matrix_ = new_matrix;
}

template <typename T>
void Matrix2D<T>::setMatrix(std::vector< std::vector<T> > & new_matrix) {
    std::size_t row, col;
    row = new_matrix.size();
    col = new_matrix[0].size();
    row_ = row;
    col_ = col;
    matrix_ = new_matrix;
}

template <typename T>
std::unique_ptr<std::vector< std::vector<T> > > Matrix2D<T>::coFactor(const std::vector< std::vector<T> > & mat,
         std::size_t row, std::size_t col, std::size_t currentSize) {

    std::unique_ptr<std::vector< std::vector<T> > > cofactor = std::make_unique<std::vector< std::vector<T> > >();
    std::size_t idx{0}, idy{0};
    
    std::size_t matSize = mat.size();
    cofactor -> resize(matSize, std::vector<T>(matSize, T()));
        
    for (std::size_t rowIndex = 0; rowIndex < currentSize; rowIndex++) {
        for (std::size_t colIndex = 0; colIndex < currentSize; colIndex++) {
            if (rowIndex != row && colIndex != col) {
                cofactor -> at(idx).at(idy++) = mat.at(rowIndex).at(colIndex);
                if (idy == currentSize - 1) {
                    idy = 0;
                    idx++;
                }
            }
        }
    }
    return std::move(cofactor);
}

template <typename T>
T Matrix2D<T>::recursiveDeterminant(const std::vector<std::vector<T> > & mat, std::size_t size) {
    T detVal{};
    if (size == 1)
        return mat.at(0).at(0);
    T sign{1};
    for (std::size_t col = 0; col < size; ++col) {
        detVal += sign*recursiveDeterminant(*coFactor(mat, 0, col, size), size - 1)*mat.at(0).at(col);
        sign = -sign;
    }
    return detVal;
}

template <typename T>
T Matrix2D<T>::det(std::string flag) {
    TIME_THIS();
    T detValue{};
    if (flag == "recursive") {
           detValue = recursiveDeterminant(matrix_, col_);
    }
    
    return detValue;
}

template <typename T>
void Matrix2D<T>::swapRows(int rowOne, int rowTwo) {
    using vecIter = typename std::vector< std::vector<T>>::iterator;
    vecIter rowOneIt = matrix_.begin();
    vecIter rowTwoIt = matrix_.begin();
    std::advance(rowOneIt, rowOne);
    std::advance(rowTwoIt, rowTwo);
    std::swap_ranges(rowOneIt -> begin(), rowOneIt -> end(), rowTwoIt -> begin());
}

template <typename T>
std::size_t Matrix2D<T>::maxInColIndex(std::size_t colIndex) {
    T maxValue = std::abs(matrix_.at(colIndex).at(colIndex));
    std::size_t maxIndex{0};
    for (std::size_t row = colIndex; row < row_; ++row) {
        if (std::abs(matrix_.at(row).at(colIndex)) >= maxValue) {
            maxValue = std::abs(matrix_.at(row).at(colIndex));
            maxIndex = row;
        }
    }
    return maxIndex;
}

template <typename T>
std::unique_ptr<std::vector<std::vector<T>>> Matrix2D<T>::getCopyPtrMatrix () {
    std::unique_ptr<std::vector< std::vector<T>>> matrixPtr = std::make_unique<std::vector< std::vector<T>>>(matrix_);
    return std::move(matrixPtr);
} 

# endif 