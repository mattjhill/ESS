#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "./GRP.hpp"

//using namespace std;

/************************/
/************************/
/*                      */
/*  STORING RESULTS     */
/*                      */
/************************/
/************************/
struct results {
public:
	int index;
    double assembleTime;
    double factorTime;
    double solveTime;
    double error;
};

/****************************/
/****************************/
/*                          */
/*  ASSEMBLE ENTIRE MATRIX  */
/*                          */
/****************************/
/****************************/
void assemble_Matrix(const int N, const int m, const Eigen::VectorXd alpha, const Eigen::VectorXd beta, const Eigen::VectorXd t, double d, Eigen::MatrixXd& A) {
	A	=	Eigen::MatrixXd(N,N);
	for (int i=0;i<N;++i) {
		A(i,i)	=	d;
		for (int j=i+1;j<N;++j) {
			A(i,j)	=	0;
			for (int p=0;p<m;++p) {
				A(i,j)+=alpha(p)*exp(-beta(p)*fabs(t(i)-t(j)));
			}
			A(j,i)	=	A(i,j);
		}
	}
}

/****************************************************/
/****************************************************/
/*                                                  */
/*  SOLVES USING USUAL PARTIAL PIVOTED LU ALGORITHM */
/*                                                  */
/****************************************************/
/****************************************************/
void usual_Algorithm(int N, int m, Eigen::VectorXd alpha, Eigen::VectorXd beta, Eigen::VectorXd t, double d, Eigen::VectorXd rhs, results& usualAlgorithm) {

    double CPS	=	CLOCKS_PER_SEC;
    clock_t start, end;

    //  The entire matrix
    Eigen::MatrixXd A;

    //  Stores the solution
    Eigen::VectorXd solution;

    //  Stores the partial pivoted LU decompostion
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_decomp;

    
    std::cout << std::endl << "Usual method..." << std::endl;
    
    //	Assembles the matrix
    start		=	clock();
    assemble_Matrix(N, m, alpha, beta, t, d, A);
    end			=	clock();
    usualAlgorithm.assembleTime =	(end-start)/CPS;
    std::cout << std::setw(30) << "Assembly time: " << std::setw(20) << 1000*usualAlgorithm.assembleTime << std::endl;
    
    
    //	Obtains the partial pivoted LU factorization
    start		=	clock();
	lu_decomp.compute(A);
    end			=	clock();
    usualAlgorithm.factorTime   =	(end-start)/CPS;
    std::cout << std::setw(30) << "Factor time: " << std::setw(20) << 1000*usualAlgorithm.factorTime << std::endl;
    
    
    //	Obtains the solution using partial pivoted LU factorization
    start		=	clock();
	solution	=	lu_decomp.solve(rhs);
    end			=	clock();
    usualAlgorithm.solveTime	=	(end-start)/CPS;
    std::cout << std::setw(30) << "Solve time: " << std::setw(20) << 1000*usualAlgorithm.solveTime << std::endl << std::endl;

    //	Error in the computed solution
    usualAlgorithm.error        =   (A*solution-rhs).cwiseAbs().maxCoeff();
    std::cout << std::setw(30) << "Maximum of ||Ax-b|| is: " << std::setw(20) << usualAlgorithm.error << std::endl << std::endl;
}

/********************************/
/********************************/
/*                              */
/*  SOLVES USING FAST ALGORITHM */
/*                              */
/********************************/
/********************************/
void fast_Algorithm(int N, int m, Eigen::VectorXd alpha, Eigen::VectorXd beta, Eigen::VectorXd t, double d, Eigen::VectorXd rhs, results& fastAlgorithm) {

    double CPS	=	CLOCKS_PER_SEC;
    clock_t start, end;

    //  Stores the solution vectors
    Eigen::VectorXd solutionFast, solex;
    
    //	Set up the ESS class
    GRP matrix(N, m, alpha, beta, t, d);
    
    std::cout << std::endl << "Fast method.." << std::endl;

    //	Assemble the matrix.
    start			=	clock();
    matrix.assemble_Extended_Matrix();
    end				=	clock();
    fastAlgorithm.assembleTime  =	(end-start)/CPS;
    std::cout << std::setw(30) << "Assembly time: " << std::setw(20) << 1000*fastAlgorithm.assembleTime << std::endl;

    //	Factorize the matrix.
    start			=	clock();
    matrix.factorize_Extended_Matrix();
    end				=	clock();
    fastAlgorithm.factorTime	=	(end-start)/CPS;
    std::cout << std::setw(30) << "Factor time: " << std::setw(20) << 1000*fastAlgorithm.factorTime << std::endl;
    
    //	Solve the linear system
    start			=	clock();
    matrix.obtain_Solution(rhs, solutionFast,solex);
    end				=	clock();
    fastAlgorithm.solveTime     =	(end-start)/CPS;
    std::cout << std::setw(30) << "Solve time: " << std::setw(20) << 1000*fastAlgorithm.solveTime << std::endl << std::endl;
    
    //	Error in the computed solution
    fastAlgorithm.error         =	matrix.obtain_Error(rhs, solex);
    std::cout << std::setw(30) << "Maximum of ||Ax-b|| is: " << std::setw(20) << fastAlgorithm.error << std::endl << std::endl;
}


/****************************/
/****************************/
/*                          */
/*  WRITE TO LATEX FIGURE   */
/*                          */
/****************************/
/****************************/
void write_To_LaTeX_Figure(const int n, const results* algorithm, std::string filename, int m) {
    std::ofstream myfile;

    std::ostringstream filenameAssembly;
    std::ostringstream filenameFactor;
    std::ostringstream filenameSolve;
    std::ostringstream filenameError;

    filenameAssembly << filename << m << "_Assembly.tex";
    filenameFactor << filename << m << "_Factor.tex";
    filenameSolve << filename << m << "_Solve.tex";
    filenameError << filename << m << "_Error.tex";

    myfile.open(filenameAssembly.str().c_str(), std::ios::out);
    for (int k=0; k<n; ++k) {
        myfile << "(" << algorithm[k].index << ", " << algorithm[k].assembleTime << ") ";
    }
    myfile.close();

    myfile.open(filenameFactor.str().c_str(), std::ios::out);
    for (int k=0; k<n; ++k) {
        myfile << "(" << algorithm[k].index << ", " << algorithm[k].factorTime << ") ";
    }
    myfile.close();
    
    myfile.open(filenameSolve.str().c_str(), std::ios::out);
    for (int k=0; k<n; ++k) {
        myfile << "(" << algorithm[k].index << ", " << algorithm[k].solveTime << ") ";
    }
    myfile.close();
    
    myfile.open(filenameError.str().c_str(), std::ios::out);
    for (int k=0; k<n; ++k) {
        myfile << "(" << algorithm[k].index << ", " << algorithm[k].error << ") ";
    }
    myfile.close();
}

/****************************************/
/****************************************/
/*                                      */
/*  SCALING WITH 'N' FOR A FIXED 'm'    */
/*                                      */
/****************************************/
/****************************************/
void scaling_With_N(int m) {
    //  Randomization
    srand(time(NULL));

    //  Get the alpha's and beta's
    Eigen::VectorXd alpha		=	Eigen::VectorXd::Ones(m)+Eigen::VectorXd::Random(m);
    Eigen::VectorXd beta		=	Eigen::VectorXd::Ones(m)+Eigen::VectorXd::Random(m);

    //  Get the diagonal entry
    double d    =	alpha.sum()+1.0;

    //  Number of N's
    int nFastAlgorithm  =   11;
    int nUsualAlgorithm =   5;

    //  Result for fast algorithm and usual algorithm
    results* fastAlgorithm  =   new results[nFastAlgorithm];
    results* usualAlgorithm =   new results[nUsualAlgorithm];

    //  System sizes
    fastAlgorithm[0].index  =   500;
    fastAlgorithm[1].index  =   1000;
    fastAlgorithm[2].index  =   2000;

    int count           =   3;

    while (count < nFastAlgorithm) {
        fastAlgorithm[count].index  =   10*fastAlgorithm[count-3].index;
        ++count;
    }

    Eigen::VectorXd rhs;
    Eigen::VectorXd time;

    for (int k=0; k<nFastAlgorithm; ++k) {
        std::cout << std::endl << "N = " << fastAlgorithm[k].index << "; m = " << m << std::endl;
        rhs     =   Eigen::VectorXd::Random(fastAlgorithm[k].index);
        time    =   10*Eigen::VectorXd::Random(fastAlgorithm[k].index);
        //  Sort the time stamps
        std::sort(time.data(), time.data()+time.size());
        if (k<nUsualAlgorithm) {
            usualAlgorithm[k].index =   fastAlgorithm[k].index;
            usual_Algorithm(fastAlgorithm[k].index, m, alpha, beta, time, d, rhs, usualAlgorithm[k]);
        }
        fast_Algorithm(fastAlgorithm[k].index, m, alpha, beta, time, d, rhs, fastAlgorithm[k]);
    }

    write_To_LaTeX_Figure(nUsualAlgorithm, usualAlgorithm, "./results/results_Usual_Algorithm_m_", m);
    write_To_LaTeX_Figure(nFastAlgorithm, fastAlgorithm, "./results/results_Fast_Algorithm_m_", m);
    
    delete[] fastAlgorithm;
    delete[] usualAlgorithm;
}

/****************************************/
/****************************************/
/*                                      */
/*  SCALING WITH 'm' FOR A FIXED 'N'    */
/*                                      */
/****************************************/
/****************************************/
void scaling_With_m(int N) {
    //  Randomization
    srand(time(NULL));

    //  Time stamps
    Eigen::VectorXd t	=	10*Eigen::VectorXd::Random(N);

    //  Sort the time stamps
    std::sort(t.data(), t.data()+t.size());

    //  Get the right hand side
    Eigen::VectorXd rhs =	Eigen::VectorXd::Random(N);

    //  Number of m's
    int mFastAlgorithm  =   20;

    //  Stores the results
    results* fastAlgorithm  =   new results[mFastAlgorithm];
    for (int k=0; k<mFastAlgorithm; ++k) {
        fastAlgorithm[k].index      =   k+1;
        std::cout << std::endl << "N = " << N << "; m = " << k+1 << std::endl;
        Eigen::VectorXd alpha		=	Eigen::VectorXd::Ones(k+1)+Eigen::VectorXd::Random(k+1);
        Eigen::VectorXd beta		=	Eigen::VectorXd::Ones(k+1)+Eigen::VectorXd::Random(k+1);

        //  Get the diagonal entry
        double d    =	alpha.sum()+1.0;
        fast_Algorithm(N, fastAlgorithm[k].index, alpha, beta, t, d, rhs, fastAlgorithm[k]);
    }

    write_To_LaTeX_Figure(mFastAlgorithm, fastAlgorithm, "./results/results_Fast_Algorithm_N_", N);
    
    delete[] fastAlgorithm;
}

/************************************************************/
/************************************************************/
/*                                                          */
/*                      THE MAIN FUNCTION                   */
/*                                                          */
/*  argv[1] = 1 implies scaling with 'N' for a fixed 'm'    */
/*  argv[1] = 2 implies scaling with 'm' for a fixed 'N'    */
/*                                                          */
/*  argv[2] = 'm' or 'N' depending on argv[1]               */
/*                                                          */
/************************************************************/
/************************************************************/
int main(int argc, char* argv[]) {
    int option  =   atoi(argv[1]);

    if (option==1) {
        //	Rank of separability is 'm'
        int m	=	atoi(argv[2]);

        scaling_With_N(m);
    }
    else if (option==2) {
        //	Number of unknowns is 'N'
        int N	=	atoi(argv[2]);

        scaling_With_m(N);
    }
}