//	GRP.hpp
//	Created by Sivaram Ambikasaran on September 2nd, 2014

#ifndef __GRP_HPP__
#define __GRP_HPP__

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//using namespace std;
//using namespace Eigen;

class GRP {
	int N;		//	Number of unknowns.
	int m;		//	Rank of the separable part.
	//	The semi-separable matrix is of the form diag(d) + triu(U*V,1) + tril((U*V)',-1).
	Eigen::VectorXd alpha;
	Eigen::VectorXd beta;
	Eigen::VectorXd t;
	Eigen::MatrixXd gamma;
	double d;	//	Diagonal entry of the matrix.
	double determinant;	//	Determinant of the extended sparse matrix.

	int M;								//	Size of the extended sparse matrix.
	int blocknnz;						//	Number of non-zeros per block.
	int nnz;							//	Total number of non-zeros.
	int nBlockSize;						//	Size of each block, will be 2m+1.

	std::vector<double> nBlockStart;                //	Starting index of each of the blocks.
	std::vector<Eigen::Triplet<double> > triplets;	//	Vector of triplets used to store the sparse matrix.
	Eigen::SparseMatrix<double> Aex;                //	The extended sparse matrix.

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > factorize;		//	Stores the factorization.

public:
	GRP(const int N, const int m, const Eigen::VectorXd alpha, const Eigen::VectorXd beta, const Eigen::VectorXd t, const double d);                        //	Constructor gets all the desired quantities.
	void assemble_Extended_Matrix();                        //	Assembles the extended sparse matrix.
	void change_Diagonal(const Eigen::VectorXd diagonal);	//	Updates the diagonal alone.
	void factorize_Extended_Matrix();                       //	Factorizes the extended sparse matrix.
	void obtain_Solution(const Eigen::VectorXd rhs, Eigen::VectorXd& solution, Eigen::VectorXd& solutionex);	//	Obtains the solution.
	double obtain_Error(const Eigen::VectorXd rhs, const Eigen::VectorXd& solex);	//	Obtain error, i.e., ||Ax-b||_{\inf}
	double obtain_Determinant();	//	Obtains the determinant of the extended sparse matrix.
};

GRP::GRP (const int N, const int m, const Eigen::VectorXd alpha, const Eigen::VectorXd beta, const Eigen::VectorXd t, const double d) {
	/************************************************/
	/*												*/
	/*	Assign all the variables inside the class.	*/
	/*												*/
	/************************************************/
	this->N		=	N;
	this->m		=	m;
	this->alpha	=	alpha;
	this->beta	=	beta;
	this->t		=	t;
	this->d		=	d;
}

void GRP::assemble_Extended_Matrix() {
	gamma		=	Eigen::MatrixXd(m,N-1);
	for (int i=0; i<m; ++i) {
		for (int j=0; j<N-1; ++j) {
			gamma(i,j)	=	exp(-beta(i)*fabs(t(j)-t(j+1)));
		}
	}
	//	Declare 2*m as a variable, since this will be used frequently.
	int twom	=	2*m;

	//	Number of non-zeros per matrix block in the extended sparse matrix. This includes the '1' on the diagonal, the two negative identity matrices above and below the diagonal, the vectors u, u', v and v'.
	blocknnz	=	6*m+1;

	//	Declare the blocksize, which is the repeating structure along the diagonal.
	nBlockSize	=	twom+1;

	//	Size of the extended sparse matrix.
	M			=	N*nBlockSize-twom;

	//	Number of non-zero entries in the matrix. The identity matrices on the supersuperdiagonals which were not included in the blocknnz has been accounted for here.
	nnz			=	(N-1)*blocknnz+(N-2)*twom+1;	//	This number is correct and has been checked with MATLAB.

	//	Obtain the starting index of each of the block.
	for (int k=0; k<N; ++k) {
		nBlockStart.push_back(k*nBlockSize);
	}

	//	Assembles block by block except the identity matrices on the supersuperdiagonal.
	for (int nBlock=0; nBlock<N-1; ++nBlock) {
		//	The starting row and column for the blocks.
		//	Assemble the diagonal first.
		triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock], nBlockStart[nBlock], d));
		for (int k=0; k<m; ++k) {
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+1,nBlockStart[nBlock],gamma(k,nBlock)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock],nBlockStart[nBlock]+k+1,gamma(k,nBlock)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+m+k+1,nBlockStart[nBlock]+twom+1,alpha(k)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+twom+1,nBlockStart[nBlock]+m+k+1,alpha(k)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+1,nBlockStart[nBlock]+k+m+1,-1.0));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+m+1,nBlockStart[nBlock]+k+1,-1.0));
		}
	}
	triplets.push_back(Eigen::Triplet<double>(M-1,M-1,d));

	//	Assebmles the supersuperdiagonal identity blocks.
	for (int nBlock=0; nBlock<N-2; ++nBlock) {
		for (int k=0; k<m; ++k) {
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+m+1,nBlockStart[nBlock]+twom+k+2,gamma(k,nBlock+1)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+twom+k+2,nBlockStart[nBlock]+k+m+1,gamma(k,nBlock+1)));
		}
	}

	//	Set the size of the extended sparse matrix.
	Aex.resize(M,M);

	//	Assemble the matrix from triplets.
	Aex.setFromTriplets(triplets.begin(), triplets.end());
}

void GRP::factorize_Extended_Matrix() {
	//	Compute the sparse LU factorization of matrix `Aex'
    factorize.compute(Aex);
}

void GRP::obtain_Solution(const Eigen::VectorXd rhs, Eigen::VectorXd& solution, Eigen::VectorXd& solutionex) {
	//	Assemble the extended right hand side `rhsex'
	Eigen::VectorXd rhsex		=	Eigen::VectorXd::Zero(M);
	for (int nBlock=0; nBlock<N; ++nBlock) {
		rhsex(nBlockStart[nBlock])	=	rhs(nBlock);
	}

	//	Obtain the solution
	solutionex	=	factorize.solve(rhsex);

	//	Desired solution vector
	solution		=	Eigen::VectorXd(N);
	for (int nBlock=0; nBlock<N; ++nBlock) {
		solution(nBlock)	=	solutionex(nBlockStart[nBlock]);
	}
}

double GRP::obtain_Error(const Eigen::VectorXd rhs, const Eigen::VectorXd& solex) {
	//	Assemble the extended right hand side `rhsex'
	Eigen::VectorXd rhsex		=	Eigen::VectorXd::Zero(M);
	for (int nBlock=0; nBlock<N; ++nBlock) {
		rhsex(nBlockStart[nBlock])	=	rhs(nBlock);
	}
	double error	=	(Aex*solex-rhsex).cwiseAbs().maxCoeff();
	return error;
}

double GRP::obtain_Determinant() {
	//	Obtain determinant
	determinant	=	factorize.logAbsDeterminant();
	return determinant;
}

#endif /* defined(__GRP_HPP__) */