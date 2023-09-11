#include "Exercise.h"

void exercise_1()
{
	int N = 10; //æÿ’Û¥Û–°

	// initialize
	std::vector<vector<double>> A(N, vector<double>(N));
	std::vector<double> b(N);
	for (int i = 0; i < N-1; i++)
	{
		A[i][i] = 2;
		A[i + 1][i] = -1;
		A[i][i + 1] = -1;
		b[i] = 15;
	}
	A[N - 1][N - 1] = 2;
	b[0] = 7;
	b[N - 1] = 14;	

	//solve
	Solve_Gauss(A, b);
}

void exercise_2_1()
{
	//settings
	int N = 10;
	
	//initialize
	std::vector<std::vector<double>> A(N, vector<double>(N));
	std::vector<double> b(N);
	for (int i = 0; i < N - 1; i++) {
		A[i][i] = 2;
		A[i][i + 1] = -1;
		A[i + 1][i] = -1;
		b[i] = 102;
	}
	A[N - 1][N - 1] = 2;
	b[0] = 202; b[N-1] = 202;
	for (int i = 0; i < N; i++) {
		//b.push_back((double)rand() / (double)RAND_MAX * 100);
//		b[i] = (double)rand() / (double)RAND_MAX * 100;
	}
	//solve
	std::cout << "2.1:" << std::endl;
	Solve_cholesky(A, b);
}

void exercise_2_2()
{
	//settings
	int N = 40;
	//initialize Hilbert matrix
	std::vector<vector<double>> A(N, std::vector<double>(N));
	std::vector<double> b(N);
	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++) {
			A[i][j] = 1.0 / (i + j + 1);
			b[i] += 1.0 / (i + j + 1);
		}
	}
	//solve
	std::cout << "Hilbert:" << std::endl;
	Solve_cholesky(A, b);
}

void exercise_3_1()
{
	//settings
	int N = 100;

	//initialize
	std::vector<std::vector<double>> A(N, vector<double>(N));
	std::vector<double> b(N);
	for (int i = 0; i < N - 1; i++) {
		A[i][i] = 100;
		A[i][i + 1] = 1;
		A[i + 1][i] = 1;
	}
	A[N - 1][N - 1] = 100;
	for (int i = 0; i < N; i++) {
		b[i] = (double)rand() / (double)RAND_MAX * 100;
	}
	//solve
	std::cout << "3.1:" << std::endl;
	Solve_Gauss(A, b);
	Solve_cholesky(A, b);
}

void exercise_3_2()
{
	//settings
	int N = 40;
	//initialize Hilbert matrix
	std::vector<vector<double>> A(N, std::vector<double>(N));
	std::vector<double> b(N);
	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++) {
			A[i][j] = 1.0 / (i + j + 1);
			b[i] += 1.0 / (i + j + 1);
		}
	}
	//solve
	std::cout << "Hilbert:" << std::endl;
	Solve_Gauss(A, b);
	Solve_cholesky(A, b);
}
