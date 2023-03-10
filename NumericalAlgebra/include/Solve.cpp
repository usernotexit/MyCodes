#pragma once
#include<iostream>
#include<vector>
#include<Windows.h>
#include<stdlib.h>
#include<time.h>

#include"Solve.h"
#include "Function.h"

clock_t start_t, end_t;
using namespace std;

void Solver::gauss_elim(){
	if (decom != 0) { return; }
	//std::cout << "Gaussian elimination ..." << std::endl;
	//start counting time
	start_t = clock();

	funs::gauss_elim(A_);
	u.clear(); for (int i = 0; i < int(A.size()); i++) { u.push_back(i); }
	v.clear(); for (int i = 0; i < int(A[0].size()); i++) { v.push_back(i); }

	//stop counting
	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
	decom = Gauss;
	lu = U;
}

void Solver::gauss_elim_full_pivoting() {
	if (decom != 0) { return; }
	//std::cout << "Gaussian full-pivoting elimination ..." << std::endl;
	start_t = clock();

	funs::gauss_elim_full_pivoting(A_,u,v);

	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
	decom = Gauss_ful;
	lu = U;
}

void Solver::gauss_elim_col_pivoting() {
	if (decom != 0) { return; }
	//std::cout << "Gaussian column-pivoting elimination ..." << std::endl;
	start_t = clock();

	funs::gauss_elim_col_pivoting(A_, u);
	v.clear(); for (int i = 0; i < int(A[0].size()); i++) { v.push_back(i); }

	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
	decom = Gauss_col;
	lu = U;
}

void Solver::cholesky_decomp() {
	if (decom != 0) { return; }
	//std::cout << "Cholesky ..." << std::endl;
	start_t = clock();

	funs::cholesky_decomp(A_);
	///u.clear(); for (int i = 0; i < int(A.size()); i++) { u.push_back(i); }
	///v.clear(); for (int i = 0; i < int(A[0].size()); i++) { v.push_back(i); }

	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
	decom = Cholesky;
	lu = U;
}

void Solver::modified_cholesky_decomp() {
	if (decom != 0) { return; }
	//std::cout << "modified cholesky ..." << std::endl;
	start_t = clock();

	funs::modified_cholesky_decomp(A_);
	funs::matrix_DLT(A_);
	u.clear(); for (int i = 0; i < int(A.size()); i++) { u.push_back(i); }
	v.clear(); for (int i = 0; i < int(A[0].size()); i++) { v.push_back(i); }

	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
	decom = modified_Cholesky;
	lu = U;
}

void Solver::solve(){
	funs::check(A_.size(), b.size());
	x.clear();
	for (unsigned int i = 0; i < A_[0].size(); i++) { x.push_back(b[i]); }

	// Cholesky
	if (decom == 4) {
		start_t = clock();
		funs::forward_subs(A_, x);		//solve Ly = b 
		funs::back_subs(A_, x);			//solve LT x = y
		end_t = clock();
		time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
		if (A.size() != 0) { funs::loss_Inf(A, x, b, loss); }
		return;
	}

	// QR
	if (decom == 6) {
		start_t = clock();
		/// TODO: Solve the n*n equations, same with Solver::solve_LS()
		int m = A_.size(), n = A_[0].size();
		if (m < n) { cout << "We can't solve this now, try new data" << endl; return; }

		// caculate QT*b
		for (int i = 0; i < A_[0].size(); i++) {
			vector<double> v;
			v.push_back(1.0);
			for (int j = i + 1; j < int(A_.size()); j++) {
				v.push_back(A_[j][i]);
			}
			funs::house_trans(x, v, beta[i], i);
		}

/*
		// get (c1, c2)T = QT*b
		vector<double> c1;
		for (int i = 0; i < n; i++) { c1.push_back(x[i]); }
*/

		// solve Rx = c1;
		funs::back_subs(A_, x);

		end_t = clock();
		time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
		if (A.size() != 0) { funs::loss_Inf(A, x, b, loss); }
		return;
	}

	// Gauss and modified_Cholesky
	switch (lu)
	{
	case 0: std::cout << "the coefficient matrix has not decomposed yet,"
			<< " please check your code and restart the program."
			<< std::endl;
		system("pause");
		exit(-1);
	case 1: start_t = clock();
		funs::vector_pb(u, x);			//calculate P*b
		funs::forward_subs(A_, x);
		funs::back_subs1(A_, x);
		funs::vector_qb(v, x);			//solve Q^-1 * x = z (x = Qz)
		end_t = clock();
		time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
		if (A.size() != 0) { funs::loss_Inf(A, x, b, loss); }
		break;
	case 2: start_t = clock();
		funs::vector_pb(u, x);			//calculate P*b
		funs::forward_subs1(A_, x);
		funs::back_subs(A_, x);
		funs::vector_qb(v, x);			//solve Q^-1 * x = z (x = Qz)
		end_t = clock();
		time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;
		if (A.size() != 0) { funs::loss_Inf(A, x, b, loss); }
		break;
	default: std::cout << "there is no matching keyword,"
		<< " please check your code and restart the program."
		<< std::endl;
		system("pause");
		exit(-2);
		break;
	}
	return;
}
void Solver::solve(vector<double> &x) {
	this->solve();
	x.clear();
	x = this->x;
}

void Solver::solve(vector<double> &x, const decom_type decom_method){

}

/*
	The following is about class 'Solver2'
*/
void Solver2::solve(int &epoch, const iter_type type ) {//= Jacobi
	x.clear();
	for (int i = 0; i < b.size(); i++) { x.push_back(b[i]); }// 初值选为b

	start_t = clock();

	int i;
	switch (type)
	{
	case Jacobi:
		for (i = 0; i < iter; i++) {
			this->updateJacobi();
			if (err < tol) { break; }
			if (i % 1000 == 0)cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
		}
		break;
	case GS:
		for (i = 0; i < iter; i++) {
			this->updateGS();
			if (err < tol) { break; }
			if (i % 1000 == 0)cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
		}
		break;
	case SOR:// factor w is required
		for (i = 0; i < iter; i++) {
			this->updateSOR();
			if (err < tol) { break; }
			if (i % 1000 == 0)cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
		}
		break;
	case CG: 
	{
		for (int i = 0; i < x.size(); i++) { x[i]=0; }//初值选为0
		vector<double> r = this->b, p = r;
		double r_2=funs::dot(r,r);
		for (i = 0; i < iter; i++) {
			this->updateCG(r, p, r_2);
			if (err < tol) { break; }
			if (i % 10 == 0) {
				cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
			}
		}
		break;
	}
	default: cout << "??" << endl;
		break;
	}

	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;	// time
	epoch = i;									// iteration times
	funs::loss_Inf(A, x, b, loss);				// loss
}

void Solver2::solve(const iter_type type){// = Jacobi
	int i;
	this->solve(i, type);
}

void Solver2::updateJacobi() {
	// B = D^-1  * (L+U)
	vector<double> y;
	funs::updateJacobi(A, b, x, y);	// update x(k+1)
	
	funs::LInfnorm(x, y, err);	//update err
	x = y;
}

void Solver2::updateGS() {
	// B = D^-1  * (L+U)
	vector<double> y;
	funs::updateGS(A, b, x, y, 1.0);

	funs::LInfnorm(x, y, err);
	x = y;
}

void Solver2::updateSOR() {
	// B = D^-1  * (L+U)
	vector<double> y;
	funs::updateGS(A, b, x, y, w);

	funs::LInfnorm(x, y, err);
	x = y;
}

void Solver2::updateCG(vector<double> &r, vector<double> &p, double &r_2) {
	/* r, p, r_2
		q = Ap
		a = r_2/(q,p)
		x += ap
		r -= aq
		tmp = (r,r)
		b = -tmp/r_2
		r_2 = tmp
		p = r+bp
	*/

	vector<double> q;
	double a, b, tmp;
	err = 0.;
	funs::times(this->A, p, q);
	a = r_2 / (funs::dot(p, q));
	//cout << a << endl;//for test
	for (unsigned int i = 0; i < x.size(); i++) {
		x[i] += a * p[i];
		if (err < fabs(a*p[i])) { err = fabs(a*p[i]); }	//update err/step
		r[i] -= a * q[i];
	}
	tmp = funs::dot(r, r);
	b = tmp / r_2;
	r_2 = tmp;
	for (unsigned int i = 0; i < p.size(); i++) {
		p[i] = r[i] + b * p[i];
	}
}

/*
	The following is about "SolverSparse" class
*/
void SolverSparse::solve(int &epoch, const iter_type type) {
	// the original x contains boundary information, so donnot cover it before iterating
	start_t = clock();

	int i;
	switch (type)
	{
	case Jacobi:
		for (i = 0; i < iter; i++) {
			this->updateJacobi();
			if (err < tol) { break; }
			if (i % 10 == 0)cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
		}
		break;
	case GS:
		for (i = 0; i < iter; i++) {
			this->updateGS();
			if (err < tol) { break; }
			if (i % 10 == 0)cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
		}
		break;
	case SOR:// factor w is required
		for (i = 0; i < iter; i++) {
			this->updateSOR();
			if (err < tol) { break; }
			if (i % 10 == 0)cout << "epoch " << i << "..." << "err:" << err << endl;// retained for test
		}
		break;
	case CG:
	{
		int m = g.size(), n = g[0].size();
		vector<vector<double>> r = f,p;//this->x
		double r_2;
		
		// folding the edge
		for (int i = 1; i <  m - 1; i++) {  // row 0, row m,column 0, column n are absorbed
			r[i][1] += x[i][0]; 
			r[i][n - 2] += x[i][n - 1];
		}
		for (int i = 1; i < n - 1; i++) {
			r[1][i] += x[0][i];
			r[m - 2][i] += x[m - 1][i];
		}
		//initialize x0 = 0
		for (int i = 1; i < m - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				x[i][j] = 0.;
			}
		}
		p = r;
		//cut down the edge: please take a look at funs::timesSparse()
		for (int i = 0; i < n; i++) {
			p[0][i] = 0; p[m - 1][i] = 0;
		}
		for (int i = 0; i < m; i++) {
			p[i][0] = 0; p[i][n - 1] = 0;
		}

		r_2 = funs::dotSparse(r, r);

		for (i = 0; i < iter; i++) {
			this->updateCG(r, p, r_2);
			if (err < tol) { break; }
		}
		break;
	}
	default: cout << "??" << endl;
		break;
	}

	end_t = clock();
	time_consuming += (double)(end_t - start_t) / CLOCKS_PER_SEC;	// time
	epoch = i;									// iteration times
	//funs::loss_Inf(A, x, b, loss);				// loss
	loss = 0.;
}

void SolverSparse::solve(const iter_type type) {
	int i;
	this->solve(i, type);
}

void SolverSparse::updateJacobi() {//update x(k+1) and err
	vector<vector<double>> y;
	funs::updateJacobiSparse(g, f, x, y);	// update x(k+1)

	funs::LInfnorm(x, y, err);	//update err
	x = y;
}

void SolverSparse::updateGS() {//update x(k+1) and err
	vector<vector<double>> y;
	funs::updateGSSparse(g, f, x, y);	// update x(k+1)

	funs::LInfnorm(x, y, err);	//update err
	x = y;
}

void SolverSparse::updateSOR() {//update x(k+1) and err
	vector<vector<double>> y;
	funs::updateSORSparse(g, f, x, y, w);	// update x(k+1)

	funs::LInfnorm(x, y, err);	//update err
	x = y;
}

void SolverSparse::updateCG(vector<vector<double>> &r, vector<vector<double>> &p, double &r_2) {
	vector<vector<double>> q; 
	double a, b, tmp;
	funs::timesSparse(g, p, q);
	a = r_2 / funs::dotSparse(p, q);
	err = 0.;
	for (unsigned int i = 1; i < r.size()-1; i++) {
		for (unsigned int j = 1; j < r[0].size() - 1; j++) {
			double t = a * p[i][j];
			r[i][j] -= a * q[i][j];
			x[i][j] += t;
			if (err < fabs(t)) err = fabs(t);	//update err
		}
	}
	tmp = funs::dotSparse(r, r);
	b = tmp / r_2;
	r_2 = tmp;
	for (unsigned int i = 1; i < r.size()-1; i++) {// 1~n-1 -> 0~n
		for (unsigned int j = 1; j < r[0].size()-1; j++) {
			p[i][j] = r[i][j] + b * p[i][j];
		}
	}
}
