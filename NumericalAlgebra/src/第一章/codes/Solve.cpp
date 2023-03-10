#pragma once
#include<iostream>
#include<vector>
#include<Windows.h>
#include<stdlib.h>
#include<time.h>
#include"Solve.h"

clock_t start_t, end_t;
using namespace std;

void Solve_Gauss(vector<vector<double>> &A0, vector<double> b0) {
	/* I would like to get the data from the files,
	but it is not convenient to maintain them */
	/*
		std::ifstream outf;
		outf.open("data/equation1.txt");
		std::infile >> N;
		...
	*/
	int N = A0.size();
	{// Gaussian  elimination
		// 初始化A和b
		auto A = A0;
		auto b = b0;

		//start counting time
		start_t = clock();
	
		std::cout << "Gaussian elimination ..." << std::endl;
		gauss_elim(A);			//decomposition A = LU
		forward_subs1(A, b);	//solve Ly = b
		back_subs(A, b);		//solve Ux = y
/*
	std::cout << "LU decomposition:" << std::endl;
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			std::cout << A[i][j] << "\t";
		}
		std::cout << "\n";
	}
*/
		//stop counting
		end_t = clock();

		std::cout << "solution:" << std::endl;
		for (int i = 0; i < N; i++) {
			std::cout << b[i] << "\t";
		}
		//show consumed time
		std::cout << "\nTime:" << (double)(end_t - start_t) / CLOCKS_PER_SEC
			<< "s" << std::endl;
		//show error
		double err = 0;
		vector_Ax(A0, b, b0, err);
		std::cout << "Error:" << err << std::endl;
	}

	{// Gaussian full-pivoting elimination
		// 初始化A和b
		auto A = A0;
		auto b = b0;
		std::vector<int> u, v;			//represents transition matrix P and Q

		//start counting time
		start_t = clock();

		std::cout << "Gaussian full-pivoting elimination ..." << std::endl;
		gauss_elim_full_pivoting(A, u, v);	//decomposition PAQ = LU
		vector_pb(u, b);			//calculate P*b
		forward_subs1(A, b);		//solve Ly = P*b (b)
		back_subs(A, b);			//solve Uz = y
		vector_qb(v, b);			//solve Q^-1 * x = z (x = Qz)

		//stop counting
		end_t = clock();

		//std::cout << "LU decomposition:" << std::endl;
		//for(int i=0; i<N; i++){
		//	for(int j=0; j<N; j++){
		//		std::cout << A[i][j] << "\t";
		//	}
		//	std::cout << "\n";
		//}
		std::cout << "solution:" << std::endl;
		for (int i = 0; i < N; i++) {
			std::cout << b[i] << "\t";
		}
		//show consumed time
		std::cout << "\nTime:" << (double)(end_t - start_t) / CLOCKS_PER_SEC
			<< "s" << std::endl;
		//show error
		double err = 0;
		vector_Ax(A0, b, b0, err);
		std::cout << "Error:" << err << std::endl;
	}

	{// Gaussian col-pivoting elimination
		//初始化A和b
		auto A = A0;
		auto b = b0;
		std::vector<int> u;			//represents transition matrix P

		//start counting time
		start_t = clock();

		std::cout << "Gaussian col-pivoting elimination ..." << std::endl;
		gauss_elim_col_pivoting(A, u);	//decomposition PA = LU
		vector_pb(u, b);			//calculate P*b
		forward_subs1(A, b);		//solve Ly = P*b (b)
		back_subs(A, b);			//solve Ux = y

		//stop counting
		end_t = clock();

		std::cout << "LU decomposition:" << std::endl;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				std::cout << A[i][j] << "\t";
			}
			std::cout << "\n";
		}
		std::cout << "solution:" << std::endl;
		for (int i = 0; i < N; i++) {
			std::cout << b[i] << "\t";
		}
		//show consumed time
		std::cout << "\nTime:" << (double)(end_t - start_t) / CLOCKS_PER_SEC
			<< "s" << std::endl;
		//show error
		double err = 0;
		vector_Ax(A0, b, b0, err);
		std::cout << "Error:" << err << std::endl;
	}
}

void Solve_cholesky(vector<vector<double>> &A0, vector<double> b0) {
	int N = A0.size();
	{//cholesky
		//initialize
		auto A = A0;
		auto b = b0;

		//solve equations
		std::cout << "cholesky ..." << std::endl;
		start_t = clock();			//time
		cholesky_decomp(A);
		forward_subs(A, b);			//solve Ly = b 
		back_subs(A, b);			//solve L^-1 x = y
		end_t = clock();

		std::cout << "LU decomposition:" << std::endl;
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				std::cout << A[i][j] << "\t";
			}
				std::cout << "\n";
		}
		//show the result
		std::cout << "solution:" << std::endl;
		for (int i = 0; i < N; i++) {
			std::cout << b[i] << "\t";
		}
		//show consumed time
		std::cout << "\nTime:" << (double)(end_t - start_t) / CLOCKS_PER_SEC
			<< "s" << std::endl;
		//show error
		double err = 0;
		vector_Ax(A0, b, b0, err);
		std::cout << "Error:" << err << std::endl;
	}

	{//modified_cholesky
		//initialize
		auto A = A0;
		auto b = b0;

		//solve equations
		std::cout << "modified cholesky ..." << std::endl;
		start_t = clock();			//time
		modified_cholesky_decomp(A);
		matrix_DLT(A);
		forward_subs1(A, b);		//solve Ly = b 
		back_subs(A, b);			//solve L^-1 x = y
		end_t = clock();

		//show the result
		std::cout << "solution:" << std::endl;
		for (int i = 0; i < N; i++) {
			std::cout << b[i] << "\t";
		}
		//show consumed time
		std::cout << "\nTime:" << (double)(end_t - start_t) / CLOCKS_PER_SEC
			<< "s" << std::endl;
		//show error
		double err = 0;
		vector_Ax(A0, b, b0, err);
		std::cout << "Error:" << err << std::endl;
	}
}

