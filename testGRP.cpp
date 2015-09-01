#include <iostream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "./GRP.hpp"

//using namespace std;

void assemble_Matrix(const int N, const int m, const Eigen::VectorXd alpha, const Eigen::VectorXd beta, const Eigen::VectorXd t, const Eigen::VectorXd d, Eigen::MatrixXd& A) {
	A	=	Eigen::MatrixXd(N,N);
	for (int i=0;i<N;++i) {
		A(i,i)	=	d(i);
		for (int j=i+1;j<N;++j) {
			A(i,j)	=	0;
			for (int p=0;p<m;++p) {
				A(i,j)+=alpha(p)*exp(-beta(p)*fabs(t(i)-t(j)));
			}
			A(j,i)	=	A(i,j);
		}
	}
}

void get_Exact_Factorization(const Eigen::MatrixXd A, Eigen::PartialPivLU<Eigen::MatrixXd>& lu_decomp) {
	//	Obtains the lu decomposition
	lu_decomp.compute(A);
}

void get_Solution(const Eigen::PartialPivLU<Eigen::MatrixXd> lu_decomp, const Eigen::VectorXd rhs, Eigen:: VectorXd& solution) {
	//	Obtains the solution
	solution	=	lu_decomp.solve(rhs);
}

int main(int argc, char* argv[]) {

	// Eigen::initParallel();
	// int n       =   atoi(argv[4]);
	// Eigen::setNbThreads(n);
	//
	//	Number of unknowns is 'N'
	int N	=	atoi(argv[1]);

	//	Rank of separability is 'm'
	int m	=	atoi(argv[2]);

	//	Do you want to check with usual solver? c=y
	char c	=	argv[3][0];

	srand(time(NULL));

	Eigen::VectorXd alpha		=	Eigen::VectorXd::Ones(m)+Eigen::VectorXd::Random(m);
	Eigen::VectorXd beta		=	Eigen::VectorXd::Ones(m)+Eigen::VectorXd::Random(m);

//  double d    =   alpha.sum();
	Eigen::VectorXd d	=	Eigen::VectorXd::Constant(N, alpha.sum()+1.0);
//	double d	=	1.0;

	Eigen::VectorXd t	=	Eigen::VectorXd::Random(N);
	std::sort(t.data(), t.data()+t.size());

	/****************************************************/
	/*													*/
	/*	Covariance kernel is of the following form		*/
	/*													*/
	/*	K(r)	=	sum_k a_k exp(-b_k |t_i-t_j|)		*/
	/*													*/
	/****************************************************/
	Eigen::VectorXd rhs=	Eigen::VectorXd::Random(N);
	Eigen::VectorXd solution;
	double start, end;

	/************************/
	/*						*/
	/*	Fast factorization	*/
	/*						*/
	/************************/


	Eigen::VectorXd solutionFast, solex;
	double assembleFastTime, factorFastTime, solveFastTime, determinantTime;

	//	Set up the ESS class.
	GRP matrix(N, m, alpha, beta, t, d);

	std::cout << std::endl << "Fast method.." << std::endl;

	//	Assemble the matrix.
	start   = omp_get_wtime();
	matrix.assemble_Extended_Matrix();
	end     = omp_get_wtime();
	assembleFastTime=	(end-start);
	std::cout << std::setw(30) << "Assembly time: " << std::setw(10) << 1000*assembleFastTime << std::endl;

	//	Factorize the matrix.
	start   = omp_get_wtime();
	matrix.factorize_Extended_Matrix();
	end     = omp_get_wtime();
	factorFastTime	=	(end-start);
	std::cout << std::setw(30) << "Factor time: " << std::setw(10) << 1000*factorFastTime << std::endl;

	//	Solve the linear system
	start   = omp_get_wtime();
	matrix.obtain_Solution(rhs, solutionFast,solex);
	end     = omp_get_wtime();
	solveFastTime	=	(end-start);
	std::cout << std::setw(30) << "Solve time: " << std::setw(10) << 1000*solveFastTime << std::endl << std::endl;

	//	Error in the computed solution
	double error 	=	matrix.obtain_Error(rhs, solex);
	std::cout << std::setw(30) << "Maximum of ||Ax-b|| is: " << std::setw(10) << error << std::endl << std::endl;

	//	Obtain the determinant of the extended sparse system
	start	=	omp_get_wtime();
	double logAbsdeterminantFast	=	matrix.obtain_Determinant();
	end		=	omp_get_wtime();
	determinantTime	=	(end-start);
	std::cout << std::setw(30) << "LogAbsDeterminant is: " << std::setw(10) << logAbsdeterminantFast << std::endl << std::endl;
	std::cout << std::setw(30) << "Determinant time: " << std::setw(10) << 1000*determinantTime << std::endl << std::endl;

	/************************************************/
	/*												*/
	/*	Direct factorization and compare solutions	*/
	/*												*/
	/************************************************/

	if (c=='y') {
		Eigen::MatrixXd A;
		Eigen::PartialPivLU<Eigen::MatrixXd> lu_decomp;

		double assembleTime, factorTime, solveTime, determinantTime;

		std::cout << std::endl << "Usual method..." << std::endl;

		//	Assembles the matrix
		start   = omp_get_wtime();
		assemble_Matrix(N, m, alpha, beta, t, d, A);
		end     = omp_get_wtime();
		assembleTime=	(end-start);
		std::cout << std::setw(30) << "Assembly time: " << std::setw(10) << 1000*assembleTime << std::endl;


		//	Obtains the factorization
		start   = omp_get_wtime();
		get_Exact_Factorization(A, lu_decomp);
		end     = omp_get_wtime();
		factorTime	=	(end-start);
		std::cout << std::setw(30) << "Factor time: " << std::setw(10) << 1000*factorTime << std::endl;


		//	Obtains the solution
		start   = omp_get_wtime();
		get_Solution(lu_decomp, rhs, solution);
		end     = omp_get_wtime();
		solveTime	=	(end-start);
		std::cout << std::setw(30) << "Solve time: " << std::setw(10) << 1000*solveTime << std::endl << std::endl;
		
		//	Error in the computed solution
		double usualError   =   (A*solution-rhs).cwiseAbs().maxCoeff();
		std::cout << std::setw(30) << "Maximum of ||Ax-b|| is: " << std::setw(10) << usualError << std::endl << std::endl;

		//	Determinant
		start	=	omp_get_wtime();

		double logAbsdeterminantExact	=	log(fabs(lu_decomp.determinant()));
		end		=	omp_get_wtime();
		determinantTime	=	(end-start);
		std::cout << std::setw(30) << "LogAbsDeterminant is: " << std::setw(10) << logAbsdeterminantExact << std::endl << std::endl;
		std::cout << std::setw(30) << "Determinant time: " << std::setw(10) << 1000*determinantTime << std::endl << std::endl;

		std::cout << std::endl << "Error in solution: "<< (solutionFast-solution).cwiseAbs().maxCoeff() << std::endl << std::endl;
		std::cout << std::endl << "Error in LogAbsDeterminant: " << fabs(logAbsdeterminantFast-logAbsdeterminantExact)/fabs(logAbsdeterminantExact) << std::endl << std::endl;
	}
}