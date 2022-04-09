# ifndef TESTER_H_
# define TESTER_H_

# include "header.h"
# include "helpers.h"

// class for testing a single function

template <typename ReturnType, typename Matrix, typename... InputArgs>

class TestFunction {

    public:
        TestFunction(InputArgs... args, std::function<ReturnType(InputArgs& ...)> testFunction):
            inputs_(std::forward_as_tuple(args...)), testFunction_{testFunction}{};
        
        template <typename... Inputs>
        TestFunction& addFixture(const std::string &, Inputs&&...);

        static std::unique_ptr<Matrix> read1DMatrix(const std::string &, auto);        
        static std::unique_ptr<Matrix> read2DMatrix(const std::string &, auto, auto);

        void run();
    
    private:
        std::tuple<InputArgs&...> inputs_;
        ReturnType outputs_;
        std::function<ReturnType(InputArgs& ...)> testFunction_;
        std::vector<Matrix> fixtures_;
};

template <typename ReturnType, typename Matrix, typename... InputArgs>
void TestFunction<ReturnType, Matrix, InputArgs...>::run() {
    ReturnType outputs = std::apply(testFunction_, inputs_);
}


template <typename ReturnType, typename Matrix, typename... InputArgs>
template <typename... Inputs>
TestFunction<ReturnType, Matrix, InputArgs...>& TestFunction<ReturnType, Matrix, InputArgs...>::addFixture(const std::string & fixtureFileName, Inputs&&... auxInputs) {
    if constexpr (Matrix::Dim == 2) {
        auto numRows = getFromTuple<0>(auxInputs...);
        auto numCols = getFromTuple<1>(auxInputs...);
        fixtures_.push_back(*read2DMatrix(fixtureFileName, numRows, numCols));
    }
    else if constexpr (Matrix::Dim == 1) {
        auto numRows = getFromTuple<0>(auxInputs...);
        fixtures_.push_back(*read1DMatrix(fixtureFileName, numRows));
    }
    return *this;
}

template <typename ReturnType, typename Matrix, typename... InputArgs>
std::unique_ptr<Matrix> TestFunction<ReturnType, Matrix, InputArgs...>::read1DMatrix(const std::string & filename, auto numRows) {

    std::ifstream input(filename);
    std::string line{};
    std::string values{};
    typedef typename Matrix::matrixType T;
    std::size_t row{0};
    std::unique_ptr<Matrix> fixture = std::make_unique<Matrix>(numRows);

    while (std::getline(input, line)) {
        std::istringstream in(line);
        while (std::getline(in, values, ' ')) {
                if constexpr (std::is_integral<T>::value) {
                    fixture -> operator()(row++) = std::stoi(values);
                }
                else if constexpr (std::is_floating_point<T>::value) {
                    fixture -> operator()(row++) = static_cast<T>(std::stod(values));
                }
        }
    }
    return std::move(fixture);
}

template <typename ReturnType, typename Matrix, typename... InputArgs>
std::unique_ptr<Matrix> TestFunction<ReturnType, Matrix, InputArgs...>::read2DMatrix(const std::string & filename, auto numRows, auto numCols) {

    std::ifstream input(filename);
    std::string line{};
    std::string values{};
    typedef typename Matrix::matrixType T;
    std::size_t row{0}, col{0};
    std::unique_ptr<Matrix> fixture = std::make_unique<Matrix>(numRows, numCols);

    while (std::getline(input, line)) {
        std::istringstream in(line);
        while (std::getline(in, values, ' ')) {
                if constexpr (std::is_integral<T>::value) {
                    fixture -> operator()(row, col++) = std::stoi(values);
                }
                else if constexpr (std::is_floating_point<T>::value) {
                    fixture -> operator()(row, col++) = std::stof(values);
                }
        }
        col = 0;
        ++row;
    }
    return std::move(fixture);
}


# endif 