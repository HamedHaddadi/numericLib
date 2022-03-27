// running basic sanity check tests on matrix class
//  until developing a test suite 
# include "../arrays/matrix2D.h"
# include "../arrays/matrix1D.h"
# include "../arrays/arrayND.h"
# include "../linearSystems/luSolvers.h"
# include "../linearSystems/choleskyDecomposition.h"

int main() {

    std::cout <<"comparing basic recursive determinant between linearized arrayND and Matrix2D ... "<<std::endl;
    std::cout <<" matrix2D >>>>"<<std::endl;
    Matrix2D<double> matOne(10,10);
    matOne.fillRandom(-4.,12.);
    double detMat = matOne.det();
    std::cout <<"det value is "<<detMat<<std::endl;
    std::cout <<" arrayND >>>>"<<std::endl;
    Array<double, 2> arrOne(10,10);
    arrOne.fillRandom(-4, 12);
    double detArr = arrOne.det();
    std::cout <<"det value is "<<detArr <<std::endl;

    std::cout <<"Testing LU decomposition "<<std::endl;
    Matrix2D<int> matDecomp(4,4);
    matDecomp.fillRandom(1.,5.);
    std::cout <<"Matrix2D is "<<std::endl;
    std::cout <<matDecomp<<std::endl;
    Matrix2D<double> U = std::get<0>(doolittleLUDecomposition(matDecomp));
    Matrix2D<double> L = std::get<1>(doolittleLUDecomposition(matDecomp));
    std::cout <<L<<std::endl;
    std::cout <<U <<std::endl;
    std::cout <<"multiply them "<<std::endl;
    std::cout <<L*U<<std::endl;

    std::cout <<"test gaussian elimination "<<std::endl;
    Matrix1D<double> b(10);
    b.fillRandom(-1.0, 20.0);
    Matrix2D<double> A(10, 10);
    A.fillRandom(-2, 14);
    Matrix1D<double> b_ans = std::get<1>(gaussianElimination(A, b));
    Matrix2D<double> A_ans = std::get<0>(gaussianElimination(A, b));
    std::cout <<"A ans is "<<std::endl;
    std::cout << A_ans<<std::endl;
    std::cout <<"b ans is "<<std::endl;
    std::cout <<b_ans <<std::endl;

    std::cout <<"testinh Cholesky decomposition on a positive definite array" <<std::endl;
    Matrix2D<double> cholTest(3,3);
    cholTest.inputFromString("4 12 -16 12 37 -43 -16 -43 98");
    std::pair<bool, Matrix2D<double> > results = choleskyDecomposition(cholTest);
    if (results.first)
        std::cout <<"cholesky exist "<<std::endl;
    std::cout <<"lower diagonal matrix is "<<std::endl;
    std::cout <<results.second;

    return 0;
}