//	ESS.hpp
//	Created by Sivaram Ambikasaran on September 2nd, 2014

#ifndef __ESS_HPP__
#define __ESS_HPP__

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//using namespace std;
//using namespace Eigen;

class ESS {
	int N;					//	Number of unknowns.
	int m;					//	Rank of the separable part.
	//	The semi-separable matrix is of the form diag(d) + triu(U*V,1) + tril((U*V)',-1).
	Eigen::MatrixXd U;
	Eigen::MatrixXd V;
	Eigen::VectorXd d;		//	Diagonal of the matrix.

	int M;								//	Size of the extended sparse matrix.
	int blocknnz;						//	Number of non-zeros per block.
	int nnz;							//	Total number of non-zeros.
	int nBlockSize;						//	Size of each block, will be 2m+1.
	std::vector<double> nBlockStart;			//	Starting index of each of the blocks.
	std::vector<Eigen::Triplet<double> > triplets;	//	Vector of triplets used to store the sparse matrix.

	Eigen::SparseMatrix<double> Aex;			//	The extended sparse matrix.
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > factorize;		//	Stores the factorization.

public:
	ESS(const int N, const int m, const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen::VectorXd diagonal);	//	Constructor gets all the desired quantities.
	void assemble_Extended_Matrix();				//	Assembles the extended sparse matrix.
	void change_Diagonal(const Eigen::VectorXd diagonal);	//	Updates the diagonal alone.
	void factorize_Extended_Matrix();				//	Factorizes the extended sparse matrix.
	void obtain_Solution(const Eigen::VectorXd rhs, Eigen::VectorXd& solution);	//	Obtains the solution.
};

ESS::ESS (const int N, const int m, const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen::VectorXd diagonal) {
	/************************************************/
	/*												*/
	/*	Assign all the variables inside the class.	*/
	/*												*/
	/************************************************/
	this->N	=	N;
	this->m	=	m;
	this->U	=	U;
	this->V	=	V;
	this->d	=	diagonal;
}

void ESS::assemble_Extended_Matrix() {
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
		triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock], nBlockStart[nBlock], d(nBlock)));
		for (int k=0; k<m; ++k) {
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+1,nBlockStart[nBlock],U(nBlock,k)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock],nBlockStart[nBlock]+k+1,U(nBlock,k)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+m+k+1,nBlockStart[nBlock]+twom+1,V(k,nBlock+1)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+twom+1,nBlockStart[nBlock]+m+k+1,V(k,nBlock+1)));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+1,nBlockStart[nBlock]+k+m+1,-1.0));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+m+1,nBlockStart[nBlock]+k+1,-1.0));
		}
	}
	triplets.push_back(Eigen::Triplet<double>(M-1,M-1,d(N-1)));

	//	Assebmles the supersuperdiagonal identity blocks.
	for (int nBlock=0; nBlock<N-2; ++nBlock) {
		for (int k=0; k<m; ++k) {
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+k+m+1,nBlockStart[nBlock]+twom+k+2,1.0));
			triplets.push_back(Eigen::Triplet<double>(nBlockStart[nBlock]+twom+k+2,nBlockStart[nBlock]+k+m+1,1.0));
		}
	}

	//	Set the size of the extended sparse matrix.
	Aex.resize(M,M);

	//	Assemble the matrix from triplets.
	Aex.setFromTriplets(triplets.begin(), triplets.end());
}

void ESS::factorize_Extended_Matrix() {
	//	Compute the sparse LU factorization of matrix `Aex'
    factorize.compute(Aex);
}

void ESS::obtain_Solution(const Eigen::VectorXd rhs, Eigen::VectorXd& solution) {
	//	Assemble the extended right hand side `rhsex'
	Eigen::VectorXd rhsex		=	Eigen::VectorXd::Zero(M);
	for (int nBlock=0; nBlock<N; ++nBlock) {
		rhsex(nBlockStart[nBlock])	=	rhs(nBlock);
	}

	//	Obtain the solution
	Eigen::VectorXd solutionex	=	factorize.solve(rhsex);

	//	Desired solution vector
	solution		=	Eigen::VectorXd(N);
	for (int nBlock=0; nBlock<N; ++nBlock) {
		solution(nBlock)	=	solutionex(nBlockStart[nBlock]);
	}
}
#endif /* defined(__ESS_HPP__) */