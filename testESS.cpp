#include <iostream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "./ESS.hpp"

//using namespace std;

void assemble_Matrix(const int N, const int m, const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen:: VectorXd d, Eigen::MatrixXd& A) {
	Eigen::MatrixXd B	=	U*V;
	//	Assembling the diagonal of the semi-separable matrix
	A	=	d.asDiagonal();
	//	Assembling the upper triangular part of the semi-separable matrix
	A.triangularView<Eigen::StrictlyUpper>()	=	B.triangularView<Eigen::StrictlyUpper>();
	//	Assembling the lower triangular part of the semi-separable matrix
	A.triangularView<Eigen::StrictlyLower>()	=	A.triangularView<Eigen::StrictlyUpper>().transpose();
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
	//	Number of unknowns is 'N'
	int N	=	atoi(argv[1]);

	//	Rank of separability is 'm'
	int m	=	atoi(argv[2]);

	//	Do you want to check with usual solver? c=y
	char c	=	argv[3][0];

	srand(time(NULL));

	Eigen::MatrixXd U	=	Eigen::MatrixXd::Random(N,m);
	Eigen::MatrixXd V	=	Eigen::MatrixXd::Random(m,N);

	Eigen::VectorXd d	=	Eigen:: VectorXd::Random(N);

	Eigen::VectorXd rhs=	Eigen:: VectorXd::Random(N);
	Eigen::VectorXd solution;

	double CPS	=	CLOCKS_PER_SEC;
	clock_t start, end;


	/****************************/
	/*							*/
	/*	Direct factorization	*/
	/*							*/
	/****************************/

	if (c=='y') {
		Eigen::MatrixXd A;

		Eigen::PartialPivLU<Eigen::MatrixXd> lu_decomp;

		double assembleTime, factorTime, solveTime;

		std::cout << std::endl << "Usual method..." << std::endl;

		//	Assembles the matrix
		start		=	clock();
		assemble_Matrix(N, m, U, V, d, A);
		end			=	clock();
		assembleTime=	(end-start)/CPS;
		std::cout << std::setw(30) << "Assembly time: " << std::setw(10) << assembleTime << std::endl;


		//	Obtains the factorization
		start		=	clock();
		get_Exact_Factorization(A, lu_decomp);
		end			=	clock();
		factorTime	=	(end-start)/CPS;
		std::cout << std::setw(30) << "Factor time: " << std::setw(10) << factorTime << std::endl;


		//	Obtains the solution
		start		=	clock();
		get_Solution(lu_decomp, rhs, solution);
		end			=	clock();
		solveTime	=	(end-start)/CPS;
		std::cout << std::setw(30) << "Solve time: " << std::setw(10) << solveTime << std::endl << std::endl;
	}


	/************************/
	/*						*/
	/*	Fast factorization	*/
	/*						*/
	/************************/


	Eigen::VectorXd solutionFast;
	double assembleFastTime, factorFastTime, solveFastTime;

	//	Set up the ESS class.
	ESS matrix(N, m, U, V, d);

	std::cout << std::endl << "Fast method.." << std::endl;

	//	Assemble the matrix.
	start			=	clock();
	matrix.assemble_Extended_Matrix();
	end				=	clock();
	assembleFastTime=	(end-start)/CPS;
	std::cout << std::setw(30) << "Assembly time: " << std::setw(10) << assembleFastTime << std::endl;

	//	Factorize the matrix.
	start			=	clock();
	matrix.factorize_Extended_Matrix();
	end				=	clock();
	factorFastTime	=	(end-start)/CPS;
	std::cout << std::setw(30) << "Factor time: " << std::setw(10) << factorFastTime << std::endl;

	//	Solve the linear system
	start			=	clock();
	matrix.obtain_Solution(rhs, solutionFast);
	end				=	clock();
	solveFastTime	=	(end-start)/CPS;
	std::cout << std::setw(30) << "Solve time: " << std::setw(10) << solveFastTime << std::endl << std::endl;

	/************************/
	/*						*/
	/*	Compare solutions	*/
	/*						*/
	/************************/

	//	Error in the computed solution
	if (c=='y') {
		std::cout << std::endl << "Error in solution: "<< (solutionFast-solution).cwiseAbs().maxCoeff() << std::endl << std::endl;		
	}
}