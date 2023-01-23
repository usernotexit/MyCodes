#pragma once
#include<iostream>
#include<vector>
#include"Function.h"

using namespace std;
/* 分解法解方程，解LS问题*/
class Solver
{
	private:
		vector<vector<double>> A;	// coefficient matrix
		vector<vector<double>> A_;	// after decompose
		vector<double> x;			// solution
		vector<int> u;				// gauss_elim, trans on rows
		vector<int> v;				// gauss_elim, trans on colunmns
		vector<double> beta;		// householder QR, beta

		decom_type decom = not_yet;	// how you decompose A
		LU_type lu = none;			// diag elements belongs to ...
		//bool symmetric = false;		// A is symmetric
		
		double time_consuming = 0.0;// 1 per second
		double loss = 0.0;			// loss: |Ax-b|Inf
		double err = 0.0;			// estimation error: max|x-x_real|Inf

	public:	
		vector<double> b;			// 
		Solver(const vector<vector<double>> &_A, const vector<double> &_b){
			A = _A; b = _b; A_ = _A; x.clear();
		} 
		Solver(const vector<vector<double>> &_A) {
			A = _A; A_ = _A; 
		}
		
		// skip decomposition step, quick equation solving
		Solver(const vector<vector<double>> &_A, const vector<double> &_b, const vector<int> &_u,
			const vector<int> &_v, const decom_type type, const LU_type type_lu) {
			A_ = _A; b = _b; x.clear();
			u = _u; for (int i = _u.size(); i < int(A.size()); i++) { u.push_back(i); }
			v = _v; for (int i = _v.size(); i < int(A.size()); i++) { v.push_back(i); }
			decom = type;
			lu = type_lu;
		}
		
		// get info or results from Solver
		void cancanA(vector<vector<double>> &_A){				// show A_
			_A.clear();
			for(unsigned int i=0; i<A_.size(); i++){ _A.push_back(A_[i]); }
		}
		void cancanA(vector<vector<double>> &_A, vector<int> &_u, vector<int> &_v) {
			_A.clear(); for (int i = 0; i < A_.size(); i++) { _A.push_back(A_[i]); }
			_u.clear(); for (int i = 0; i < u.size(); i++) { _u.push_back(u[i]); }
			_v.clear(); for (int i = 0; i < v.size(); i++) { _v.push_back(v[i]); }
		}
		void showshow_result(vector<double> &x) { x = this->x; }
		void showshow_result(vector<double> &x, double &time, double &loss) {
			x = this->x;
			time = this->time_consuming;
			loss = this->loss;
		}
		void showshow_result(vector<double> &x, double &time, double &loss, double &err) {
			x = this->x;
			time = this->time_consuming;
			loss = this->loss;
			err = this->err;
		}
		void cancanA(vector<vector<double>> &_A, vector<double> &_v) { _A = A_; _v = beta; }

		// Decomposition Methods
		void gauss_elim();
		void gauss_elim_full_pivoting();
		void gauss_elim_col_pivoting();
		void cholesky_decomp();
		void modified_cholesky_decomp();
		void QR();
		
		// solve equations with A decomposed into A_
		void solve(vector<double> &x);
		void solve();

		// solve equations, you can choose the method
		void solve(vector<double> &x, const decom_type decom_method);
		// analyse the conditional number of A
		void con_num(double &con_num);
		// analyse estimation error: max|x-x_real|Inf
		void error();

		// solve LS problem
		void solve_LS();
		void solve_LS(double &min);
		void solve_LS(vector<double> &x, double &min);
}; 

/*
	迭代法解方程
*/
class Solver2
{
	private:
		vector<vector<double>> A;	// coefficient matrix
		vector<double> b;			// 
		vector<double> x;			// solution
		int iter = 1e5;				// maximum 迭代次数
		double tol = 1e-6;			// tolerance 停机条件
		double w = 1.2;				// relaxation factor in SOR

		double time_consuming = 0.0;// 1 per second
		double loss = 0.0;			// loss: |Ax-b|Inf
		double err = 0.0;			// error (or step): max|x(k)-x(k+1)|Inf

	public:
		Solver2(const vector<vector<double>> &_A, const vector<double> &_b) {
			A = _A; b = _b; x.clear();
		}
		Solver2(const vector<vector<double>> &_A, const vector<double> &_b, double _w) {
			A = _A; b = _b; x.clear(); w = _w;
		}
		Solver2(const vector<vector<double>> &_A, const vector<double> &_b, const double _tol , const int _iter){
			A = _A; b = _b; x.clear(); 
			tol = _tol; iter = _iter;
		}
		Solver2(const vector<vector<double>> &_A, const vector<double> &_b, const double _tol, const int _iter, double _w) {
			A = _A; b = _b; x.clear();
			tol = _tol; iter = _iter;
			w = _w;
		}

		// get info or results from Solver
		void show_result(vector<double> &x) { x = this->x; }
		void show_result(vector<double> &x, double &time, double &loss) {
			x = this->x;
			time = this->time_consuming;
			loss = this->loss;
		}
		void show_result(vector<double> &x, double &time, double &loss, double &err) {
			x = this->x;
			time = this->time_consuming;
			loss = this->loss;
			err = this->err;
		}

		// solve the equations by iterative methods
		void solve(int &iter, const iter_type type = Jacobi);
		void solve(const iter_type type = Jacobi);

		// update x each step
		void updateJacobi();
		void updateGS();
		void updateSOR();
		void updateCG(vector<double> &r, vector<double> &p, double &r_2);
};

/*
	针对离散化二维拉普拉斯方程结构的求解器——具有特定结构的稀疏矩阵
	boundary conditions are required
	-(detla)u + g(x,y) * u = f(x,y)
*/
class SolverSparse
{
private:
	vector<vector<double>> g;	// coefficient matrix
	vector<vector<double>> f;	// 
	vector<vector<double>> x;	// solution
	int iter = 1e5;				// maximum 迭代次数
	double tol = 1e-7;			// tolerance
	double w = 1.2;				// relaxation factor in SOR

	double time_consuming = 0.0;// 1 per second
	double loss = 0.0;			// loss: |Ax-b|Inf
	double err = 0.0;			// error (or step): max|x(k)-x(k+1)|Inf
public:
	SolverSparse(vector<vector<double>> &_g, vector<vector<double>> &_f, vector<vector<double>> &_x, int _iter=1e5) {
		g = _g; f = _f; iter = _iter; x = _x;
	}
	SolverSparse(vector<vector<double>> &_g, vector<vector<double>> &_f, vector<vector<double>> &_x, double _w) {
		g = _g; f = _f; w = _w; x = _x;
	}
	SolverSparse(vector<vector<double>> &_g, vector<vector<double>> &_f, vector<vector<double>> &_x, double _w, double _tol,int _iter=1e5) {
		g = _g; f = _f; w = _w; tol = _tol; iter = _iter; x = _x;
	}

	void show_result(vector<vector<double>> &x) { x = this->x; }
	void show_result(vector<vector<double>> &x, double &time, double &loss) {
		x = this->x;
		time = this->time_consuming;
		loss = this->loss;
	}
	void show_result(vector<vector<double>> &x, double &time, double &loss, double &err) {
		x = this->x;
		time = this->time_consuming;
		loss = this->loss;
		err = this->err;
	}

	// solve the equations by iterative methods
	void solve(int &iter, const iter_type type = Jacobi);
	void solve(const iter_type type = Jacobi);

	// update x each step
	void updateJacobi();
	void updateGS();
	void updateSOR();
	void updateCG(vector<vector<double>> &r, vector<vector<double>> &p, double &r_2);

};