/*
	This code is about QR decomposition, solving linear equations and LS problems by this method.
*/
#pragma once

#include<iostream>
#include<vector>
#include <cmath>
#include<time.h>

#include "Function.h"
#include "Solve.h"
clock_t start_QR, end_QR;

void Solver::QR() {
	if (decom != 0) { return; }
	start_QR = clock();

	funs::QR_decom(A_, beta);

	end_QR = clock();
	time_consuming += (double)(end_QR - start_QR) / CLOCKS_PER_SEC;
	decom = decom_type(6);
}

void Solver::solve_LS() {
	start_QR = clock();

	if (decom == 0) { this->QR(); }
	if (decom != 6) { cout << "please do QR decomposition" << endl; return; }

	funs::check(A_.size(), b.size());
	int m = A_.size(), n = A[0].size();
	if (m < n) { cout << "We can't solve this now, try new data" << endl; return; }
	x.clear();
	for (int i = 0; i < b.size(); i++) { x.push_back(b[i]); }


	// caculate QT*b
	for (int i = 0; i < A_[0].size(); i++) {
		vector<double> v;
		v.push_back(1.0);
		for (int j = i + 1; j < int(A_.size()); j++) {
			v.push_back(A_[j][i]);
		}
		funs::house_trans(x, v, beta[i], i);
	}

	// get (c1, c2)T = QT*b
	for (int i = m; i > n; i--) { x.pop_back(); }

	// solve Rx = c1;
	funs::back_subs(A_, x);

	end_QR = clock();
	time_consuming += (double)(end_QR - start_QR) / CLOCKS_PER_SEC;
}

void Solver::solve_LS(double &min) {
	this->solve_LS();
	funs::loss_2(A, x, b, min);
}

void Solver::solve_LS(vector<double> &x, double &min){
	this->solve_LS();
	x = this->x;
	funs::loss_2(A, this->x, b, min);
}

