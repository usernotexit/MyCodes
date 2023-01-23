/*
	Algorithm 2.5.1: caculate L1 norm of matrix A^-1, furthermore, the conditional number of A
*/
#pragma once

#include<iostream>
#include<vector>
#include <cmath>

#include "Function.h"
#include "Solve.h"

using namespace std;

// calculate LInf norm of A
void LInfnorm(vector<vector<double>> &A, double &norm) {
	norm = 0.0;
	for (unsigned int i = 0; i < A.size(); i++) {
		double tmp = 0;
		for (auto k : A[i]) {
			tmp += fabs(k);
		}
		if (tmp > norm) {
			norm = tmp;
		}
	}
}

//calculate LInf norm of b
void LInfnorm(const vector<double> b, double &norm, int &j) {
	norm = 0.0; j = 0; int j1 = 0;
	for (auto k : b) {
		if (fabs(k) > norm) { norm = fabs(k); j = j1; }
		j1++;
	}
}

//calculate LInf norm of b
void LInfnorm(const vector<double> b, double &norm) {
	norm = 0.0;
	for (auto k : b) {
		if (fabs(k) > norm) { norm = fabs(k); }
	}
}

// 最优化 估计矩阵A^-1的无穷范数 和对应的向量x, A has been decomposed into LU form
// cholesky decomposition is not avaliable now
void LInv_Infnorm(vector<vector<double>> &A, double &norm, vector<double> &x,
	vector<int> &A_u, vector<int> &A_v, const decom_type type){
	
	if(! A.size()){ cout << "error!" << endl; return;}
	
	x.clear(); for(unsigned int i=0; i<A.size(); i++){ x.push_back(1.0/A.size()); }
	bool is_min = false;
	vector<double> w;
	vector<double> z;
	vector<double> v(A.size());
	vector<vector<double>> AT(A.size(), vector<double>(A.size()));
	for (unsigned int i = 0; i < A.size(); i++) {		// transform
		for (unsigned int j = 0; j < A.size(); j++) {
			AT[i][j] = A[j][i];
		}
	}

	while(!is_min){
		// w = (A^-1)T x; v = sign(w)
		Solver S1(AT, x, A_v, A_u, type, L);	/// in need of improvement
		S1.solve(w);
		for (unsigned int i = 0; i < A.size(); i++) {
			v[i] = double(w[i] > 0) - double(w[i] < 0);
		}
		
		// z = (A^-1)*v
		Solver S2(A, v, A_u, A_v, type, U);		/// in need of improvement
		S2.solve(z);

		// ||z|| vs zT*x
		double normInf_z = 0.0, zx=0.0; int k;
		LInfnorm(z, normInf_z, k);
		for(unsigned int i=0; i<z.size(); i++){ zx += x[i]*z[i]; }
		
		if (normInf_z <= fabs(zx)) { norm = normInf_z; is_min = true; }	// areadly minimized
		else{
			x.clear();
			for(unsigned int i=0; i<A.size(); i++){x.push_back(0);}
			x[k] = 1;
		}
	}
}

void Solver::con_num(double & con_num) {
	if (decom == 0) { this->gauss_elim(); }
	//if(decom !=)
	double _A, _A_Inv;
	vector<double> t;
	LInfnorm(A, _A);
	LInv_Infnorm(A_, _A_Inv, t, u, v, decom);
	con_num = _A * _A_Inv;
}

void Solver::error() {
	if (x.size() == 0) {
		if (decom == 0) {
			std::cout << "automatically decompose ..." << std::endl;
			this->gauss_elim_col_pivoting();
		}
		this->solve();
	}
	double norm_r, norm_b, con_num;
	this->con_num(con_num);
	norm_r = loss;
	LInfnorm(b, norm_b);
	this->err = con_num * norm_r / norm_b;
}