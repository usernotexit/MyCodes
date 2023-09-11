#include <iostream>
#include<vector>
#include<time.h>

#include "Err.h"

using namespace std;

int main()
{
	//exerise 1
	cout << "exercise 1:" << endl;
	for (int N = 5; N <= 30; N++) {
		double con_num;
		vector<vector<double>> A(N,vector<double>(N));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i][j] = 1.0 / (i + j + 1);
			}
		}
		Solver s(A);
		s.gauss_elim();
		//s.gauss_elim_col_pivoting();
		s.con_num(con_num);		// gauss elimiation
		cout << "N = " << N << ": " << con_num << endl;
	}


	//exercise 2
	cout << "exercise 2:" << endl;
	for (int N = 5; N <= 30; N++) {
		double time,err,Loss,loss,true_err,norm_b;
		vector<vector<double>> A(N, vector<double>(N));
		vector<double> b(N); vector<double> x(N),x_;
		for (int i = 0; i < N; i++) {
			A[i][i] = 1.0;
			for (int j = 0; j < i; j++) {
				A[i][j] = -1.0;
				A[j][i] = 0.0;
			}
		}
		for (int i = 0; i < N; i++) {
			A[i][N-1] = 1.0;
		}
		for (int i = 0; i < N; i++) {
			x[i] = double(rand())* 1000/7;
		}
		
		// caculate b = Ax
		for (int i = 0; i < N; i++) {
			b[i] = 0;
			for (int j = 0; j < N; j++) {
				b[i] += A[i][j] * x[j];
			}
		}

		// Solve by gauss_col_privoting method, the result is x_
		Solver s(A,b);
		s.gauss_elim_col_pivoting();
		s.error();
		s.showshow_result(x_, time, Loss, err);

		// caculate true_err
		for (int i = 0; i < N; i++) { x_[i] -= x[i]; }
		LInfnorm(b, norm_b);
		LInfnorm(x_, loss);
		true_err = loss / norm_b;
		
		cout << "N = " << N << ": "
			<< "time: " << time * 1000.0 << "ms\t"
			<< "loss: " << Loss << "\t"
			<< "err: " << err << '\t'
			<< "true_err: " << true_err << endl;

	}
}

