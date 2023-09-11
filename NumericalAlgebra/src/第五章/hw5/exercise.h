#pragma once
#include<iomanip>
#include "Solve.h"

using namespace std;

void exe1_CG(const int n, vector<vector<double>> g, vector<vector<double>> f, vector<vector<double>> phi,
	double tol);
void exe1_SOR(const int n, vector<vector<double>> g, vector<vector<double>> f, vector<vector<double>> phi,
	double tol, double w);
void exe2_Hilbert(const int n);
void exe3_solve(Solver2 s, iter_type type);
void print_matrix(vector<vector<double>> x);
void print_vector(vector<double> x);

void exe1_CG(const int n, vector<vector<double>> g, vector<vector<double>> f, vector<vector<double>> phi,
	double tol) {
	cout << "CG..." << endl;
	int iter;
	double time, loss;// loss is meaningless here
	vector<vector<double>> x;
	SolverSparse s(g, f, phi, 1, tol);
	s.solve(iter, CG);

	s.show_result(x, time, loss);
	print_matrix(x);// 
	cout << "time: " << time << "s" << endl;
	cout << "iteration numbers: " << iter << endl;
}

void exe1_SOR(const int n, vector<vector<double>> g, vector<vector<double>> f, vector<vector<double>> phi,
	double tol, double w) {
	cout << "SOR..." << endl;
	cout << "w = " << w << endl;
	int iter;
	double time, loss;// loss is meaningless here
	vector<vector<double>> x;
	SolverSparse s(g, f, phi, w, tol);
	s.solve(iter, SOR);

	s.show_result(x, time, loss);
	//print_matrix(x);//option
	cout << "time: " << time << "s" << endl;
	cout << "iteration numbers: " << iter << endl;
}

void exe2_Hilbert(const int n) {
	int iter; double time, loss;
	vector<vector<double>> A(n, vector<double>(n));
	vector<double> b(n), x;
	for (int i = 0; i < n; i++) {
		b[i] = 0;
		for (int j = 0; j < n; j++) {
			A[i][j] = 1.0 / (i + j + 1);
			b[i] += A[i][j] / 3;
		}
	}
	Solver2 s(A, b, 1e-7, 1000, 1.7);//tol=1e-7
	s.solve(iter, SOR);
	s.show_result(x, time, loss);
	print_vector(x);//option
	cout << "time: " << time << "s" << endl;
	cout << "iteration numbers: " << iter << endl;

}

void print_matrix(vector<vector<double>> x) {
	cout.precision(4);					// 设置实数显示四位小数
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < x[0].size(); j++) {
			cout << x[i][j] << '\t';
		}
		cout << endl;
	}
	cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
}

void print_vector(vector<double> x) {
	cout.precision(4);					// 设置实数显示四位小数
	for (int i = 0; i < x.size(); i++) {
		cout << x[i] << '\t';
	}
	cout << endl;
	cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
}

void exe3_solve(Solver2 s, iter_type type) {
	double time, loss, err;
	int iter;
	vector<double> c;
	switch (type)
	{
	case CG:cout << "CG..." << endl; break;
	case Jacobi:cout << "Jacobi..." << endl; break;
	case GS:cout << "GS..." << endl; break;
	default:cout << "err" << endl;
		break;
	}

	s.solve(iter, type);
	s.show_result(c, time, loss);
	print_vector(c);
	cout << "time: " << time*1000.0 << "ms" << endl;
	cout << "iteration numbers: " << iter << endl;
}