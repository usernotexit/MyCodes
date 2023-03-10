#pragma once
#include<iostream>
#include<iomanip>

#include "Solve.h"

using namespace std;

/*
	epsilon = 
	a = 1/2
	n = 100
	h = 1/N
*/
void exercise_1(int N, double a, double epsilon, double tol) {
	// initialize
	std::vector<vector<double>> A(N-1, vector<double>(N-1)), B;
	std::vector<double> b(N-1), x;
	double time, loss, h=1.0/N;
	int iter;

	for (int i = 0; i < N - 2; i++)
	{
		A[i][i] = -(2 * epsilon + h);
		A[i + 1][i] = epsilon;
		A[i][i + 1] = epsilon + h;
		b[i] = a*h*h;
	}
	A[N - 2][N - 2] = -(2 * epsilon + h);
	b[0] = a * h*h;
	b[N - 2] = a * h*h - (epsilon + h)*1.0;

	{//Jacobi
		cout << "Jacobi..." << endl;
		Solver2 S(A, b, tol, 100000);
		S.solve(iter, Jacobi);
		S.show_result(x, time, loss);		// loss is meaningless here

		loss = 0;
		for (int i = 0; i < x.size(); i++) {
			double tmp = (1 - a)*(1 - exp(-(i+1)*h/(double)epsilon)) / (1-exp(-1.0/epsilon)) 
				+ a*(i+1)*h;
			if (loss < fabs(x[i] - tmp)) loss = fabs(x[i] - tmp);
			//cout << fabs(x[i] - tmp)<<endl;
		}
		/**/
		cout << "result: " << endl;
		cout.precision(4);					// 设置实数显示四位小数
		for (int i = 0; i < x.size(); i++) { cout << x[i] << "\t"; }
		cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
		
		cout << "\niterate times: " << iter << endl;
		cout << "time_consuming: " << time << "s" << endl;
		cout << "bias: " << loss << endl;
	}

	{//G-S
		cout << "G-S..." << endl;
		Solver2 S(A, b, tol, 100000);
		S.solve(iter, GS);
		S.show_result(x, time, loss);		// loss is meaningless here

		loss = 0;
		for (int i = 0; i < x.size(); i++) {
			double tmp = (1 - a)*(1 - exp(-(i + 1)*h / (double)epsilon)) / (1 - exp(-1.0 / epsilon))
				+ a * (i + 1)*h;
			if (loss < fabs(x[i] - tmp)) loss = fabs(x[i] - tmp);
		}
		/**/
		cout << "result: " << endl;
		cout.precision(4);					// 设置实数显示四位小数
		for (int i = 0; i < x.size(); i++) { cout << x[i] << "\t"; }
		cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
		
		cout << "\niterate times: " << iter << endl;
		cout << "time_consuming: " << time << "s" << endl;
		cout << "bias: " << loss << endl;
	}

	for (int k = 0; k < 3; k++)
	{//SOR
		cout << "SOR..." << " w=" << 1.1 + 0.1*k << endl;
		Solver2 S(A, b, tol, 100000, 1.1+0.1*k);
		S.solve(iter, SOR);
		S.show_result(x, time, loss);		// loss is meaningless here

		loss = 0;
		for (int i = 0; i < x.size(); i++) {
			double tmp = (1 - a)*(1 - exp(-(i + 1)*h / (double)epsilon)) / (1 - exp(-1.0 / epsilon))
				+ a * (i + 1)*h;
			if (loss < fabs(x[i] - tmp)) loss = fabs(x[i] - tmp);
		}
		/**/
		cout << "result: " << endl;
		cout.precision(4);					// 设置实数显示四位小数
		for (int i = 0; i < x.size(); i++) { cout << x[i] << "\t"; }
		cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
		
		cout << "\niterate times: " << iter << endl;
		cout << "time_consuming: " << time << "s" << endl;
		cout << "bias: " << loss << endl;
	}
}

/*
	g = e^xy, f = x+y 
	tol = 1e-7 
	edge = 1 
*/
void exercise_2(int N) {
	std::vector<vector<double>> g((N+1), vector<double>(N+1)), f((N+1), vector<double>(N+1));
	std::vector<vector<double>> u((N + 1), vector<double>(N + 1));
	const double tol = 1e-7;
	double time, loss;
	int iter;

	for (int i = 0; i < N+1; i++)
	{
		for (int j = 0; j < N+1; j++) {
			g[i][j] = exp((i / (double)N)*(j / (double)N))/N/N;
			f[i][j] = ((i + j) / (double)N)/N/N;
		}
	}
	//egde 
	for (int i = 0; i < N+1; i++) {
		for (int j = 0; j < N+1; j++) {
			u[i][j] = 1.;
		}
	}
	
	// Jacobi & G-S
	for (int k = 0; k < 2; k++)
	{
		std::vector<vector<double>> x;
		if(k==0) cout << "Jacobi..." << endl;
		if (k == 1) cout << "G-S..." << endl;
		SolverSparse S(g, f, u);
		S.solve(iter, (iter_type)k);
		S.show_result(x, time, loss);		// loss is meaningless here

		/**/
		cout << "result: " << endl;
		cout.precision(4);					// 设置实数显示四位小数
		for (int i = 0; i < x.size(); i++) {
			for (int j = 0; j < x[0].size(); j++) {
				cout << x[i][j] << "\t"; 
			}
			cout << endl;
		}
		
		cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
		cout << "\niterate times: " << iter << endl;
		cout << "time_consuming: " << time << "s" << endl;
		// cout << "bias: " << loss << endl;
	}

	// SOR
	for (int k = 0; k < 10; k++)
	{
		std::vector<vector<double>> x;
		cout << "SOR..." << " w=" << 1.05 + 0.05*k << endl;
		double w = 1.05 + 0.05*k;
		SolverSparse S(g, f, u, w);

		S.solve(iter, SOR);
		S.show_result(x, time, loss);		// loss is meaningless here

		/**/
		cout << "result: " << endl;
		cout.precision(4);					// 设置实数显示四位小数
		for (int i = 0; i < x.size(); i++) {
			for (int j = 0; j < x[0].size(); j++) {
				cout << x[i][j] << "\t";
			}
			cout << endl;
		}
		
		cout << resetiosflags(ios::left | ios::fixed);  //清除状态左对齐和定点格式
		cout << "\niterate times: " << iter << endl;
		cout << "time_consuming: " << time << "s" << endl;
	}
}