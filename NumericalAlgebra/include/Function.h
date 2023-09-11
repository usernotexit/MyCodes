#pragma once
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

enum decom_type { not_yet, Gauss, Gauss_ful, Gauss_col, Cholesky, modified_Cholesky, QR };
enum LU_type { none, L, U };
enum iter_type { Jacobi, GS, SOR, CG };

namespace funs {
	// 矩阵运算
	void LInfnorm(vector<vector<double>> &A, double &norm);// 矩阵A的无穷范数
	void LInfnorm(vector<vector<double>> &A, vector<vector<double>> &B, double &norm);// 矩阵A-B的无穷范数
	void LInfnorm(const vector<double> b, double &norm, int &j);//向量b的无穷范数
	void LInfnorm(const vector<double> b, double &norm);//向量b的无穷范数
	void LInfnorm(const vector<double> a, const vector<double> b, double &norm);// a-b 的无穷范数
	void loss_Inf(vector<vector<double>> A, vector<double> x, vector<double> b, double &detla);//计算偏差 Ax-b
	void loss_2(vector<vector<double>> A, vector<double> x, vector<double> b, double &loss);// 计算偏差 Ax-b
	void dot(vector<double> x, vector<double> y, double &result);//向量x和y的内积
	void times(vector<vector<double>> A, vector<double> x, vector<double> &y);//计算矩阵乘法A*x
	double dot(vector<double> x, vector<double> y);//计算向量内积，直接返回结果
	void minus(vector<double> x, vector<double> y, vector<double> &z);// 向量x-y=z
	void dot(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C); // 矩阵乘法C=A B
	void minus(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C); // 矩阵减法C=A-B
	void transpose(vector<vector<double>> &A); //矩阵转置
	// householder
	void householder(vector<double> x, vector<double> &v, double &b);//计算Householder分解，归一化
	void house_trans_l(vector<vector<double>> &A, vector<double> v, double b, int i_start);//计算矩阵A的householder变换（左乘）
	void house_trans_r(vector<vector<double>> &A, vector<double> v, double b, int j_start);//计算矩阵A的householder变换（右乘）
	// givens
	void givens_trans_r(vector<vector<double>> &A, const int p, const int q, const double c, const double s);// 右乘Givens
	void givens_trans_l(vector<vector<double>> &A, const int p, const int q, const double c, const double s);// 左乘Givens

	// 三角阵 解方程
	void forward_subs(vector<vector<double>>& L, vector<double>& b);//前代法
	void forward_subs1(vector<vector<double>>& L, vector<double>& b);//对角元为1的前代法
	void back_subs(vector<vector<double>>& U, vector<double>& y);//回代法
	void back_subs1(vector<vector<double>>& U, vector<double>& y);//对角元为1的回代法
	void vector_pb(vector<int>&u, vector<double>&b);//计算向量P*b【可选】
	void vector_qb(vector<int>& v, vector<double>& b);//计算向量Q*b【可选】
	void matrix_DLT(vector<vector<double>>& A);//计算矩阵D*L^T【可选】
	// 分解法
	void gauss_elim(vector<vector<double>>& A);//Gauss消去法
	void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//全主元Gauss消去法
	void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法
	void cholesky_decomp(vector<vector<double>>& A);//对称正定阵标准Cholesky分解
	void modified_cholesky_decomp(vector<vector<double>>& A);//改进的平方根法
	// QR分解
	void QR_decom(vector<vector<double>> &A, vector<double> &v);	//QR分解矩阵
	void house_trans(vector<vector<double>> &A, vector<double> v, const double b, const int i, const int j);//计算对矩阵A的householder变换
	void house_trans(vector<double> &b, vector<double> v, const double beta, const int j_start);//计算对b的househoer变换，注意具体情境的要求
	// 迭代法
	void updateJacobi(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y);// 计算Jacobi迭代的一次更新
	void updateGS(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y, double w);// 计算G-S和SOR迭代的一次更新
	void updateJacobiSparse(const vector<vector<double>> g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y);//计算Jacobi迭代法的一次更新
	void updateGSSparse(const vector<vector<double>> g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y);//计算G-S迭代法的一次更新
	void updateSORSparse(const vector<vector<double>> g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y, const double w);//计算SOR迭代法的一次更新
	void updateCGSparse(const vector<vector<double>>g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y );//计算共轭梯度法的一次更新
	double dot(vector<vector<double>> x, vector<vector<double>> y);//计算内积x*y，x，y拉直
	void timesSparse(vector<vector<double>> g, vector<vector<double>> p, vector<vector<double>> &q);//计算矩阵乘法A*p，A由g决定
	double dotSparse(vector<vector<double>> x, vector<vector<double>> y);//计算内积x*y，x，y拉直 // specially for Sparse

	void check(int sizeL, int sizeb);	//检查矩阵行数与方程个数是否相符 
	void check(double a);				//检查对角元是否为 0 

}


namespace out {
	// output系列函数专门用于输出结果，方便展示
	void output(vector<double> a, int pre=6);

	void output(vector<vector<double>> A, int pre=6);
}
