// 实对称阵求特征值和特征向量的方法
// 过关Jacobi方法
// 二分法 + 反幂法
#pragma once
#include<iostream>
#include<vector>
#include<random>
#include"Solve.h"
using namespace std;

// 过关Jacobi方法，其中Givens变换可以单独取出
class Threshold_Jacobi
{
public:
	double thres_min;
	Threshold_Jacobi(double thres_min_= 1e-15) { this->thres_min = thres_min_; }

	void Jacobi(vector<vector<double>> A, vector<double> &v, vector<vector<double>> &V) {
		const int n = A.size();
		if (n == 0 || n != A[0].size()) { printf("本函数只处理非空方阵，请检查输入矩阵！\n"); return; }

		v.clear(); V.clear();
		// 将V初始化为对角阵
		Init(V, n);

		double thres;			//初始关值,最终关值
		_calculateE(A, thres);
		
		for (; thres>thres_min; thres/=n) {	// 过关十轮
			///printf("%f", thres);
			bool stepdown = false;
			while (!stepdown) {
				stepdown = true;
				for (int p = 1; p < n; p++) {		// 逐层扫描
					for (int q = 0; q < p; q++) {
						if (fabs(A[p][q]) > thres) { 
							double c, s;
							_givens(A, p, q, c, s);
							givens_trans2(A, p, q, c, s); //对A和V做givens变换
							givens_trans_r(V, p, q, c, s);
							///cout << "A_test" << endl;
							///out::output(A);
							stepdown = false; 
						}
					}
				}
			}
		}
		///out::output(A);
		for (int i = 0; i < n; i++) {	// 获取特征值
			v.push_back(A[i][i]);
		}
	}

	// 对矩阵A做givens相似变换(p,q,c,s)
	void givens_trans2(vector<vector<double>> &A, const int p, const int q, const double c, const double s) {
		const int n = A.size();
		if (p >= n || q >= n) { printf("err：givens变换索引超出矩阵范围！"); return; }

		// 右乘
		for (int i = 0; i < n; i++) {
			double tmp_ip = A[i][p];
			double tmp_iq = A[i][q];
			A[i][p] = c * tmp_ip - s * tmp_iq;
			A[i][q] = s * tmp_ip + c * tmp_iq;
		}
		// 左乘
		for (int i = 0; i < n; i++) {
			double tmp_pi = A[p][i];
			double tmp_qi = A[q][i];
			A[p][i] = c * tmp_pi - s * tmp_qi;
			A[q][i] = s * tmp_pi + c * tmp_qi;
		}
	}

	// 右乘Givens
	void givens_trans_r(vector<vector<double>> &A, const int p, const int q, const double c, const double s) {
		const int n = A.size();
		if (p >= n || q >= n) { printf("err：givens变换索引超出矩阵范围！"); return; }
		// 右乘
		for (int i = 0; i < n; i++) {
			double tmp_ip = A[i][p];
			double tmp_iq = A[i][q];
			A[i][p] = c * tmp_ip - s * tmp_iq;
			A[i][q] = s * tmp_ip + c * tmp_iq;
		}
	}

private:
	void _givens(vector<vector<double>> A, const int p, const int q, double &c, double &s) {
		double tau = (A[q][q] - A[p][p]) / 2 / A[p][q];
		double tan;
		tan = tau >= 0 ? (1 / (tau + sqrt(1 + tau * tau))) : (-1 / (-tau + sqrt(1 + tau * tau)));
		c = sqrt(1 / (1 + tan * tan));
		s = tan*c;
		///givens_trans2(A, p, q, c, s);
	}

	void _calculateE(vector<vector<double>> A, double &E) {
		const int n = A.size();
		E = 0.;
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				E += A[i][j] * A[i][j];
			}
		}
		E = sqrt(E);
	}

	void Init(vector<vector<double>> &V, int n) {
		for (int i = 0; i < n; i++) {
			vector<double> vi(n);
			for (int j = 0; j < n; j++) {
				if (i == j) { vi[j] = 1.; }
				else { vi[j] = 0.; }
			}
			V.push_back(vi);
		}
	}
};

// 二分法
class divided
{
public:
	// 初始化，pre指目标精度
	divided(double pre=1e-10) {
		if (pre > 0) { this->u = pre; }
		else { this->u = 1e-10; }
	}

	// 获取矩阵的第k个特征值及其对应的（某个）特征向量
	void eigensystem(vector<vector<double>> A, const int k, double &lamb, vector<double> &v) {
		this->eigenvalue(A, k, lamb);
		this->Inv_pow_method(A, lamb, v);
	}

	// 获取矩阵的第k和特征值
	void eigenvalue(vector<vector<double>> A, const int k, double &lamb) {
		double up, low, tmp;
		int s_tmp;
		funs::LInfnorm(A, up); low = -up;
		
		// 存储三对角阵的对角线和次对角线
		_get_diag(A);

		// 二分法迭代
		while (up - low > u) {
			tmp = (up + low) / 2;
			this->_get_s(tmp, s_tmp);
			s_tmp >= k ? (up = tmp) : (low = tmp);
			///printf("%d\n",s_tmp);
		}
		lamb = (up + low) / 2;
		///alpha.clear(); beta_2.clear();
	}

	// 反幂法求特征向量（请确保lamb是A的特征值）
	void Inv_pow_method(vector<vector<double>> A, const double lamb, vector<double> &v) {
		const int n = A.size();
		double lamb_;
		// u也可以充当位移精度
		for (int i = 0; i < n; i++) { A[i][i] -= lamb ; }
		
		Solver S(A);
		S.gauss_elim_col_pivoting();
		///S.cancanA(A);
		///cout << "A\n";
		///out::output(A);
/**/
		// 随机数
		vector<double> z;
		for (int k = 0; k < n; k++) { z.push_back((rand()+round(100*lamb))/1000.); }
		///out::output(z);

		// 固定若干次反幂法
		for (int epoch = 0; epoch < 10; epoch++) {
			double norm = 0;
			S.b.clear();
			S.b = z;
			S.solve(z);
			///out::output(z);
			for (int iz = 0; iz < n; iz++) { if (fabs(norm) < fabs(z[iz])) norm = z[iz]; }
			if (fabs(norm) < 1) printf("Warning: the norm is too tiny.\n");
			///printf("%f\t", norm);
			for (int iz = 0; iz < n; iz++) { z[iz] = z[iz] / norm; if(fabs(z[iz])<1e-16)z[iz]=0.;}
		}
		v = z;
	}

private:
	vector<double> alpha;		// 存储三对角阵的对角线
	vector<double> beta_2;		// 次对角线的平方，存起来方便计算
	double u;

	void _get_diag(vector<vector<double>> A) {
		const int n = A.size();
		alpha.clear(); beta_2.clear();
		for (int i = 0; i < n - 1; i++) {
			alpha.push_back(A[i][i]);
			beta_2.push_back(A[i + 1][i] * A[i + 1][i]);
		}
		if (n > 0) { alpha.push_back(A[n - 1][n - 1]); }
	}

	void _get_s(double x, int &s) {
		if (alpha.size() < 1) { printf("please check the input!\n"); s = 0; return; }
		double q;
		const int n = alpha.size();

		s = 0;
		q = alpha[0] - x;
		for (int i = 0; i < n; i++) {
			if (q < 0) { s++; }
			if (i < n-1) {
				if (q == 0) { q = sqrt(beta_2[i])*u; }
				q = alpha[i + 1] - x - beta_2[i] / q;
			}
		}
	}
};