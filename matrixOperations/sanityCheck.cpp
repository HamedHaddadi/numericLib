// running basic sanity check tests on matrix class
//  until developing a test suite 
# include "../arrays/matrix2D.h"
# include "../arrays/matrix1D.h"
# include "../arrays/arrayND.h"
# include "../linearSystems/luSolvers.h"
# include "../linearSystems/choleskyDecomposition.h"
# include "../linearSystems/jacobi.h"
# include "../utils/tester.h" 


int addThis (int & a, int & b) {
    int c{};
    c = a + b;
    return c;
}

std::pair<int, int> justReturn(int & a, int & b) {
    return std::make_pair(a, b);
}

int main() {

/*
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
    gaussianElimination(A, b);
    std::cout <<"A ans is "<<std::endl;
    std::cout << A<<std::endl;
    std::cout <<"b ans is "<<std::endl;
    std::cout <<b <<std::endl;
    */

    std::cout <<"testing Cholesky decomposition on a positive definite array" <<std::endl;
    Matrix2D<double> cholTest(3,3);
    cholTest.inputFromString("4 12 -16 12 37 -43 -16 -43 98");
    std::pair<bool, Matrix2D<double> > results = choleskyDecomposition(cholTest);
    if (results.first)
        std::cout <<"cholesky exist "<<std::endl;
    std::cout <<"lower diagonal matrix is "<<std::endl;
    std::cout <<results.second;

    std::cout <<"testing tester class "<<std::endl;
    TestFunction<std::pair<bool, Matrix2D<double>>, Matrix2D<double>, Matrix2D<double>> testThis(cholTest, &choleskyDecomposition<double>);
    testThis.addFixture("../testFixtures/cholesky3by3.dat", 3, 3);

    std::cout <<" testing jacobi method " <<std::endl;
    Matrix2D<double> jacTestA(3,3);
    Matrix1D<int> jacTestb(3);
    Matrix1D<double> jacTestRes(3);
    jacTestA.inputFromString("5 -2 3 -3 9 1 2 -1 -7");
    std::cout <<jacTestA<<std::endl;
    jacTestb.inputFromString("-1 2 3");
    std::cout << jacTestb <<std::endl;
    jacTestRes = jacobi(jacTestA, jacTestb, 1000);

    std::cout <<jacTestRes <<std::endl;


    return 0;
}