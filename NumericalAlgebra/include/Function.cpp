#include "Function.h"
#include <cmath>

using namespace std;

namespace funs {

	void forward_subs(vector<vector<double>>& L, vector<double>& b)
	{
		// 前代法解方程 Lx = b，将结果存在变量b内，L是下三角阵，b是列向量 
		check(L.size(), b.size());
		for (unsigned int k = 0; k < b.size(); k++) {
			// 从前往后逐个求解 
			b[k] = b[k] / L[k][k];
			for (unsigned int j = k + 1; j < b.size(); j++) {
				b[j] = b[j] - b[k] * L[j][k];		// b[k] has been replaced by x[k]
			}
		}
	}

	void forward_subs1(vector<vector<double>>& L, vector<double>& b)
	{
		// 前代法解方程 Lx = b，将结果存在变量b内，其中L是对角元为1的下三角阵，b是列向量，
		// 此时 L 可能不存储1，因此使用的数据范围在对角线以下（不含对角线） 
		check(L.size(), b.size());
		for (unsigned int k = 0; k < b.size(); k++) {
			for (unsigned int j = k + 1; j < b.size(); j++) {
				b[j] = b[j] - b[k] * L[j][k];
			}
		}
	}

	void back_subs(vector<vector<double>>& U, vector<double>& y)
	{
		// 回代法解方程 Ux = y，结果存在向量b内，U是上三角阵，b是列向量
		// check(U.size(), y.size());
		for (int k = int(y.size()) - 1; k > -1; k--) {
			// 从后往前逐个求解 
			y[k] = y[k] / U[k][k];
			for (int j = k - 1; j > -1; j--) {
				y[j] = y[j] - y[k] * U[j][k];		// b[k] has been replaced by x[k]
			}
		}
	}

	void back_subs1(vector<vector<double>>& U, vector<double>& y)
	{
		// 回代法解方程 Ux = y，结果存在向量b内，U是对角元为1的上三角阵，b是列向量
		// 此时 U 可能不存储1，因此使用的数据范围在对角线以下（不含对角线） 
		check(U.size(), y.size());
		for (int k = y.size() - 1; k > -1; k--) {
			y[k] = y[k] / U[k][k];
			for (int j = k - 1; j > -1; j--) {
				y[j] = y[j] - y[k] * U[j][k];		// b[k] has been replaced by x[k]
			}
		}
	}

	void gauss_elim(vector<vector<double>>& A)
	{
		// 高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
		// 参考教材上的算法
		for (unsigned int k = 0; k < A.size(); k++) {
			check(A[k][k]);
			// 计算L的第(k+1)列 
			for (unsigned int j = k + 1; j < A.size(); j++) {
				A[j][k] = A[j][k] / A[k][k];
			}
			// 让剩下部分适应前面的变换，自动计算U的第(k+1)行 
			for (unsigned int i = k + 1; i < A.size(); i++) {
				for (unsigned int j = k + 1; j < A.size(); j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
	}

	void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
	{
		// 全主元高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
		// 参考教材上的算法，PAQ = LU  
		u.clear();
		v.clear();
		for (unsigned int k = 0; k < A.size(); k++) {

			// 选取主元 
			int p = k, q = k;		//记录主元位置 
			for (unsigned int i = k; i < A.size(); i++) {
				for (unsigned int j = k; j < A.size(); j++) {
					if (fabs(A[i][j]) > fabs(A[p][q])) {
						p = i;
						q = j;
					}
				}
			}
			for (unsigned int i = 0; i < A.size(); i++) {	// exchange row k with row p 
				double swamp;
				swamp = A[k][i];
				A[k][i] = A[p][i];
				A[p][i] = swamp;
			}
			for (unsigned int j = 0; j < A.size(); j++) {	//exchange column k with column q
				double swamp;
				swamp = A[j][k];
				A[j][k] = A[j][q];
				A[j][q] = swamp;
			}
			u.push_back(p);
			v.push_back(q);

			// 接下来和普通的高斯消元法一样 
			check(A[k][k]);
			for (unsigned int j = k + 1; j < A.size(); j++) {
				A[j][k] = A[j][k] / A[k][k];
			}
			for (unsigned int i = k + 1; i < A.size(); i++) {
				for (unsigned int j = k + 1; j < A.size(); j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
	}

	void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
	{
		// 列主元高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
		// 参考教材上的算法，PA = LU 
		u.clear();
		for (unsigned int k = 0; k < A.size(); k++) {

			// 选取主元 
			int p = k;		//记录主元位置 
			for (unsigned int i = k; i < A.size(); i++) {
				if (fabs(A[i][k]) > fabs(A[p][k])) {
					p = i;
				}
			}
			for (unsigned int i = 0; i < A.size(); i++) {	// exchange row k with row p 
				double swamp;
				swamp = A[k][i];
				A[k][i] = A[p][i];
				A[p][i] = swamp;
			}

			u.push_back(p);

			// 接下来和普通的高斯消元法一样 
			check(A[k][k]);
			for (unsigned int j = k + 1; j < A.size(); j++) {
				A[j][k] = A[j][k] / A[k][k];
			}
			for (unsigned int i = k + 1; i < A.size(); i++) {
				for (unsigned int j = k + 1; j < A.size(); j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
	}

	void vector_pb(vector<int>& u, vector<double>& b)
	{
		// 行变换 P*b 
		// u represents matrix P

		// check(u, b.size());
		for (unsigned int k = 0; k < b.size(); k++) {
			double swap;
			swap = b[k];
			b[k] = b[u[k]];
			b[u[k]] = swap;
		}
	}

	void vector_qb(vector<int>& v, vector<double>& b)
	{
		// 列变换 Q*z
		// v represents matrix Q

		// check(v, b.size());

		for (unsigned int k = b.size(); k >= 1; k--) {	// the order is just reversed
			double swap;
			swap = b[k - 1];
			b[k - 1] = b[v[k - 1]];
			b[v[k - 1]] = swap;

		}
	}

	void cholesky_decomp(vector<vector<double>>& A)
	{
		// 对正定对称阵用cholesky算法分解，把结果储存在矩阵的下三角部分 
		for (unsigned int k = 0; k < A.size(); k++) {
			A[k][k] = sqrt(A[k][k]);
			for (unsigned int j = k + 1; j < A.size(); j++) {
				A[j][k] = A[j][k] / A[k][k];
			}
			for (unsigned int j = k + 1; j < A.size(); j++) {
				for (unsigned int i = j; i < A.size(); i++) {
					A[i][j] = A[i][j] - A[i][k] * A[j][k];
				}
			}
		}

		// 补上上三角部分
		for (unsigned int i = 0; i < A.size(); i++) {
			for (unsigned int j = i + 1; j < A.size(); j++) {
				A[i][j] = A[j][i];
			}
		}
	}

	void modified_cholesky_decomp(vector<vector<double>>& A)
	{
		// 对正定对称阵用改进的cholesky算法分解，把结果储存在矩阵的下三角部分 
		for (unsigned int j = 0; j < A.size(); j++) {
			std::vector<double> v(A.size());
			for (unsigned int i = 0; i < j; i++) {
				v[i] = A[j][i] * A[i][i];
				A[j][j] = A[j][j] - A[j][i] * v[i];
			}
			for (unsigned int i = j + 1; i < A.size(); i++) {
				for (unsigned int k = 0; k < j; k++) {
					A[i][j] = A[i][j] - A[i][k] * v[k];
				}
				A[i][j] = A[i][j] / A[j][j];
			}
		}

		// 补上上三角部分
		for (unsigned int i = 0; i < A.size(); i++) {
			for (unsigned int j = i + 1; j < A.size(); j++) {
				A[i][j] = A[j][i];
			}
		}
	}

	void matrix_DLT(vector<vector<double>>& A)
	{
		for (unsigned int i = 0; i < A.size(); i++) {
			for (unsigned int j = i + 1; j < A.size(); j++) {
				A[i][j] = A[i][j] * A[i][i];
			}
		}
	}

	void loss_Inf(vector<vector<double>> A, vector<double> x, vector<double> b, double &detla) {
		//计算偏差 |Ax-b|Inf，结果存在detla中
		detla = 0.0;
		check(A.size(), x.size());
		for (unsigned int i = 0; i < A.size(); i++) {
			for (unsigned int j = 0; j < A[i].size(); j++) {
				b[i] -= A[i][j] * x[j];
			}
			if (detla < fabs(b[i])) { detla = fabs(b[i]); }
		}
	}

	void loss_2(vector<vector<double>> A, vector<double> x, vector<double> b, double &loss) {
		// 计算偏差 |Ax-b|2 ，结果存在loss中
		loss = 0.0;
		check(A.size(), x.size());
		for (unsigned int i = 0; i < A.size(); i++) {
			double tmpi = b[i];
			for (unsigned int j = 0; j < A[i].size(); j++) {
				tmpi -= A[i][j] * x[j];
			}
			loss += tmpi * tmpi;
		}
		loss = sqrt(loss);
	}

	/*
		void solve_loss(vector<vector<double>> A, vector<double> x, vector<double> b, vector<double> &loss) {
			//计算偏差 Ax-b，结果存在loss中
			check(A.size(), x.size());
			loss.clear();
			for (unsigned int i = 0; i < A.size(); i++) {
				loss.push_back(b[i]);
				for (unsigned int j = 0; j < A[i].size(); j++) {
					loss[i] -= A[i][j] * x[j];
				}
			}
		}
	*/

	void LInfnorm(vector<vector<double>> &A, double &norm) {
		// calculate LInf norm of A
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
	void LInfnorm(vector<vector<double>> &A, vector<vector<double>> &B, double &norm) {
		check(A.size());
		check(A.size(), B.size());
		norm = fabs(A[0][0] - B[0][0]);
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < A[0].size(); j++) {
				if (norm < fabs(A[i][j] - B[i][j]))norm = fabs(A[i][j] - B[i][j]);
			}
		}
	}
	void LInfnorm(const vector<double> b, double &norm, int &j) {
		//calculate LInf norm of b
		norm = 0.0; j = 0; int j1 = 0;
		for (auto k : b) {
			if (fabs(k) > norm) { norm = fabs(k); j = j1; }
			j1++;
		}
	}
	void LInfnorm(const vector<double> b, double &norm) {
		//calculate LInf norm of b
		norm = 0.0;
		for (auto k : b) {
			if (fabs(k) > norm) { norm = fabs(k); }
		}
	}
	void LInfnorm(const vector<double> a, const vector<double> b, double &norm) {
		vector<double> c;
		for (int i = 0; i < a.size(); i++) { c.push_back(a[i] - b[i]); }
		funs::LInfnorm(c, norm);
	}

	void dot(vector<double> x, vector<double> y, double &result) {
		check(x.size(), y.size());
		result = 0;
		for (unsigned int i = 0; i < x.size(); i++) { result += x[i] * y[i]; }
	}
	void times(vector<vector<double>> A, vector<double> x, vector<double> &y) {
		y.clear();
		if ((A.size() == 0) || (A[0].size() == 0)) return;
		check(A[0].size(), x.size());
		for (unsigned int i = 0; i < A.size(); i++) {
			double tmp = 0;
			for (unsigned int j = 0; j < A[0].size(); j++) {
				tmp += A[i][j] * x[j];
			}
			y.push_back(tmp);
		}
	}
	double dot(vector<double> x, vector<double> y) { double result; dot(x, y, result); return result; }
	void minus(vector<double> x, vector<double> y, vector<double> &z) {
		check(x.size(), y.size());
		z.clear();
		for (unsigned int i=0; i < x.size(); i++) { z.push_back(x[i] - y[i]); }
	}
	void minus(vector<double> x, vector<double> y, vector<double> &z) {
		check(x.size(), y.size());
		z.clear();
		for (unsigned int i = 0; i < x.size(); i++) { z.push_back(x[i] - y[i]); }
	}
	void dot(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C) {
		const int m = A.size(), n = A[0].size(), l = B[0].size();
		C.clear();
		vector<double> c(l);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < l; j++) {
				c[j] = 0.;
				for (int k = 0; k < n; k++) {
					c[j] += A[i][k] * B[k][j];
				}
			}
			C.push_back(c);
		}
	}
	void minus(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C) {
		const int m = A.size(), n = A[0].size();
		C.clear(); C = A;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] -= B[i][j];
			}
		}
	}
	void transpose(vector<vector<double>> &A) {
		const int m = A.size(), n = A[0].size();
		vector<vector<double>> B(n, vector<double>(m));
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				B[j][i] = A[i][j];
			}
		}
		A = B;
	}

	void householder(vector<double> x, vector<double> &v, double &b) {
		// 输入向量x，计算householder变换Q使得Qx只有第一个分量不为0，其中Q用v和b表示，详见教材算法
		// Q = I - b v*vT
		// v的第一分量固定为1，不保留
		int n = x.size(); if (!n) { cout << "err:householder" << endl; return; }
		v.clear();
		double normInf, halfnorm, a, v0_tmp;
		LInfnorm(x, normInf);
		if(normInf <1e-15){
			for (int i = 1; i < n; i++) { v.push_back(x[i]); } 
			b = 0;
			return;
		}
		v0_tmp = x[0] / normInf;
		x[0] = v0_tmp;
		for (int i = 1; i < n; i++) { x[i] = x[i] / normInf; v.push_back(x[i]); }
		halfnorm = 0;
		for (int i = 1; i < x.size(); i++) { halfnorm += x[i] * x[i]; }

		if (halfnorm <= 1e-16) { b = 0; }
		else
		{
			a = sqrt(x[0] * x[0] + halfnorm);
			if (x[0] <= 0) { v0_tmp = x[0] - a; }
			else { v0_tmp = -halfnorm / (x[0] + a); }
			b = 2 * v0_tmp * v0_tmp / (halfnorm + v0_tmp * v0_tmp);
			for (unsigned int i = 1; i < n; i++) { v[i - 1] = v[i - 1] / v0_tmp; }	// n-1 elements in vector v
		}
	}
	void house_trans_l(vector<vector<double>>& A, vector<double> v, double b, int i_start)
	{
		const int m = A.size(), n = A[0].size();
		for (int j = 0; j < n; j++) {
			// w = ajT*v
			double w = A[i_start][j];
			for (int i = i_start + 1; i <= i_start + v.size(); i++) { w += A[i][j] * v[i - i_start - 1]; }
			// bj = aj - b*w*v
			A[i_start][j] -= b * w;
			for (int i = i_start + 1; i <= i_start + v.size(); i++) { A[i][j] -= b * w * v[i - i_start - 1]; }
		}
	}
	void house_trans_r(vector<vector<double>>& A, vector<double> v, double b, int j_start)
	{
		const int m = A.size(), n = A[0].size();
		for (int i = 0; i < m; i++) {
			// w = aiT*v
			double w = A[i][j_start];
			for (int j = j_start + 1; j <= j_start + v.size(); j++) { w += A[i][j] * v[j - j_start - 1]; }
			// bi = ai - b*w*v
			A[i][j_start] -= b * w;
			for (int j = j_start + 1; j <= j_start + v.size(); j++) { A[i][j] -= b * w * v[j - j_start - 1]; }
		}
	}
	void house_trans(vector<vector<double>> &A, vector<double> v, const double b, const int i_start, const int j_start) {
		// 计算A经过一轮householder变换后的矩阵R
		// 注意v是m-1-j_start维的，省略了第一个分量1；w是n-i_start维的
		vector<double> w;
		check(A.size()); check(A[0].size());
		int m = A.size(), n = A[0].size();

		// w = b*AT*v
		for (int i = i_start; i < n; i++) {
			double tmp = A[j_start][i];
			for (int j = j_start + 1; j < m; j++) {
				tmp += A[j][i] * v[j - j_start - 1];
			}
			w.push_back(b*tmp);
		}

		// B = A - v*wT
		for (int i = i_start; i < n; i++) {
			A[j_start][i] -= w[i - i_start];
			for (int j = j_start + 1; j < m; j++) {
				A[j][i] -= w[i - i_start] * v[j - j_start - 1];
			}
		}

	}
	void house_trans(vector<double> &b, vector<double> v, const double beta, const int j_start) {
		// 计算b经过一轮householder变换后的结果
		// 注意v是m-j_start维的
		double k;
		int m = b.size();
		if (m < j_start + 1) { cout << "something wrong with funs::house_trans." << endl; return; }

		// k = vT*b
		k = 0.0; for (int j = j_start; j < m; j++) { k += b[j] * v[j - j_start]; }
		// b' = b - beta k *v
		for (int j = j_start; j < m; j++) { b[j] -= beta * k*v[j - j_start]; }
	}

	/* givens */
	void givens_trans_r(vector<vector<double>> &A, const int p, const int q, const double c, const double s) {
		const int m = A.size(), n=A[0].size();
		///if (p >= n || q >= n) { printf("err：givens变换索引超出矩阵范围！"); return; }
		if (fabs(c*c + s * s - 1) > 1e-5) { printf("Warning!\n"); }
		// 右乘
		for (int i = 0; i < m; i++) {
			double tmp_ip = A[i][p];
			double tmp_iq = A[i][q];
			A[i][p] = c * tmp_ip - s * tmp_iq;
			A[i][q] = s * tmp_ip + c * tmp_iq;
		}
	}
	void givens_trans_l(vector<vector<double>> &A, const int p, const int q, const double c, const double s) {
		const int m = A.size(), n = A[0].size();
		///if (p >= m || q >= m) { printf("err：givens变换索引超出矩阵范围！"); return; }
		if (fabs(c*c + s * s - 1) > 1e-5) { printf("Warning!\n"); }
		// 左乘
		for (int i = 0; i < n; i++) {
			double tmp_pi = A[p][i];
			double tmp_qi = A[q][i];
			A[p][i] = c * tmp_pi + s * tmp_qi;
			A[q][i] = -s * tmp_pi + c * tmp_qi;
		}
	}

	void QR_decom(vector<vector<double>> &A, vector<double> &d) {
		// 对行数为m，列数为n的矩阵A进行QR分解，结果储存在A的上三角部分和下半部分，每一次变换的beta存在数组d中
		// 注意使用前面两条函数前。要先进行适应性调整（例如增加起始位置(i,j)以便在原矩阵A上截取子矩阵）
		check(A.size()); check(A[0].size());
		d.clear();
		int m = A.size(), n = A[0].size();
		for (int j = 0; j < n; j++) {
			double b;
			vector<double> x, v;
			for (int i = j; i < m; i++) { x.push_back(A[i][j]); }

			householder(x, v, b);
			//for (int i = j + 1; i < m; i++) { cout << v[i - 1] << '\t'; }cout <<"next:"<< endl;
			house_trans(A, v, b, j, j);
			d.push_back(b);
			for (int i = j + 1; i < m; i++) {
				A[i][j] = v[i - j - 1];
			}

		}
	}

	void check(int sizeL, int sizeb) {
		//检查矩阵行数与方程个数是否相符 
		if (sizeL != sizeb) {
			printf("Error! Please check the size of matrix.\n");
			exit(-1);
		}
	}
	void check(double a) {
		//检查对角元是否为 0 
		if (fabs(a) <= 0) {//1e-16
			printf("Cannot solve the equation. Reason: the diagonal is too close to 0.\n");
			exit(1);
		}
	}

	/* iteration steps */
	void updateJacobi(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y) {
		// y = ((L+U)*x + b) * D^-1
		y.clear();

		for (int i = 0; i < A.size(); i++) {
			double tmp = 0;
			for (int j = 0; j < A[0].size(); j++) {
				if (j == i)continue;
				else tmp -= A[i][j] * x[j];
			}
			tmp += b[i];
			funs::check(A[i][i]);
			y.push_back(tmp / A[i][i]);
		}
		//cout << y.size() << endl;
	}

	void updateGS(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y, double w) {
		// D * y = (1-w)x + w * (L*y + U*x + b)
		y.clear();
		for (int i = 0; i < x.size(); i++) { y.push_back(x[i]); }

		for (int i = 0; i < A.size(); i++) {
			double tmp = 0;
			for (int j = 0; j < A[0].size(); j++) {
				if (j == i) continue;
				else tmp -= A[i][j] * y[j];
			}
			tmp += b[i];
			y[i] = (1 - w) * y[i] + w * tmp / A[i][i];
		}
	}

	void updateJacobiSparse(const vector<vector<double>> g, const vector<vector<double>> f,
			const vector<vector<double>> x, vector<vector<double>> &y) {
		vector<vector<double>> tmp = x;// (x.size(), vector<double>(x[0].size()));

		for (int i = 1; i < x.size()-1; i++) {
			for (int j = 1; j < x[0].size()-1; j++) {
				tmp[i][j] = (f[i][j] + x[i - 1][j] + x[i][j - 1] + x[i + 1][j] + x[i][j + 1]) / (4 + g[i][j]);
			}
		}
		y = tmp;
	}
	void updateGSSparse(const vector<vector<double>> g, const vector<vector<double>> f,
		vector<vector<double>> x, vector<vector<double>> &y) {
		vector<vector<double>> tmp = x;// (x.size(), vector<double>(x[0].size()));

		for (int i = 1; i < x.size() - 1; i++) {
			for (int j = 1; j < x[0].size() - 1; j++) {
				tmp[i][j] = (f[i][j] + tmp[i - 1][j] + tmp[i][j - 1] + x[i + 1][j] + x[i][j + 1]) / (4 + g[i][j]);
			}
		}
		y = tmp;
	}
	void updateSORSparse(const vector<vector<double>> g, const vector<vector<double>> f,
		vector<vector<double>> x, vector<vector<double>> &y, const double w) {
		vector<vector<double>> tmp = x;// (x.size(), vector<double>(x[0].size()));

		for (int i = 1; i < x.size() - 1; i++) {
			for (int j = 1; j < x[0].size() - 1; j++) {
				tmp[i][j] = (1-w)*tmp[i][j] + w * (f[i][j] + tmp[i - 1][j] + tmp[i][j - 1] + x[i + 1][j] + x[i][j + 1]) / (4 + g[i][j]);
			}
		}
		y = tmp;
	}
	void updateCGSparse(const vector<vector<double>>g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y ) {

	}
	double dot(vector<vector<double>> x, vector<vector<double>> y) {
		funs::check(x.size(), y.size());
		funs::check(x[0].size(), y[0].size());
		double result = 0;
		for (unsigned int i = 0; i < x.size(); i++) {
			for (unsigned int j = 0; j < x[0].size(); j++) {
				result += x[i][j] * y[i][j];
			}
		}
		return result;
	}

	void timesSparse(vector<vector<double>> g, vector<vector<double>> p, vector<vector<double>> &q) {
		int m = g.size(), n = g[0].size();
		vector<vector<double>> tmp(m,vector<double>(n));
		for (unsigned int i = 1; i < m-1; i++) {
			for (unsigned int j = 1; j < n-1; j++) {//若不能把边界算入，则先把边界裁去（赋值0），此处可以不改
				tmp[i][j] = (4 + g[i][j])*p[i][j] - p[i - 1][j] - p[i][j - 1] - p[i + 1][j] - p[i][j + 1];
			}
		}
		q = tmp;
	}
	double dotSparse(vector<vector<double>> x, vector<vector<double>> y) {//计算内积x*y，x，y拉直 // specially for Sparse
		funs::check(x.size(), y.size());
		funs::check(x[0].size(), y[0].size());
		double result = 0;
		for (unsigned int i = 1; i < x.size()-1; i++) {
			for (unsigned int j = 1; j < x[0].size()-1; j++) {
				result += x[i][j] * y[i][j];
			}
		}
		return result;
	}
}

namespace out {
	// output系列函数专门用于输出结果，方便展示
	void output(vector<double> a, int pre) {
		const int n = a.size();
		cout << setprecision(pre);		/// 设置显示精度
		for (int i = 0; i < n; i++) {
			cout << a[i] << '\t';
			if (i % 10 == 9) { cout << '\n'; }	// 十数换行
		}
		cout << endl;
		cout << setprecision(6);		//调回默认值
	}

	void output(vector<vector<double>> A, int pre) {
		// 输出矩阵，或者一组数组
		const int n = A.size();
		cout << setprecision(pre);		/// 设置显示精度
		for (int i = 0; i < n; i++) {
			int m = A[i].size();
			for (int j = 0; j < m; j++) {
				// 较小的数输出为0，看起来更加美观
				// (fabs(A[i][j])>1e-8?(cout << A[i][j] << '\t'):(cout<< 0 <<'\t'));
				cout << A[i][j] << '\t';
			}
			cout << '\n';
		}
		cout << endl;
		cout << setprecision(6);		//调回默认值
	}
}