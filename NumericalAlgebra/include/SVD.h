/* SVD 分解算法 */
/*
	算法流程：
	二对角化
	收敛性判定及QR迭代
		后r阶：已收敛到对角阵；前l阶：待处理；中间n-l-r阶子阵H22是迭代的对象
		QR迭代：若H22有对角元为0，则调用_adapt()将对应行全部变成0
		  否则调用_qr()执行迭代
*/
#pragma once
#include<vector>
#include"Function.h"

using namespace std;
class SVD
{
public:
	SVD() { this->u = 1e-8; }
	SVD(double pre = 1e-8) { pre > 0 ? (this->u = pre) : (this->u = 1e-8); }

	//对矩阵A做SVD迭代，A=P*S*Q，结果存储在三个二维数组中
	void decom(vector<vector<double>> A, vector<vector<double>> &P, 
		vector<vector<double>> &S, vector<vector<double>> &Q);

private:
	double u;
	// 二对角化 P*S*Q = B
	void _hessenberg(vector<vector<double>> &P,
		vector<vector<double>> &S, vector<vector<double>> &Q);
	// 对S22做一步qr迭代，并将变换矩阵累乘到P和Q上
	void _qr(vector<vector<double>> &P, vector<vector<double>> &S,
		vector<vector<double>> &Q, const int r, const int l);
	// 将较小的非对角元置零
	void _set_zero(vector<vector<double>> &S);
	// 检查S到对角阵的收敛进度
	void _check_con(vector<vector<double>> &S, int &m, int &l);
	// 将P、Q初始化为单位阵
	void _init(vector<vector<double>> &V, int n) {
		for (int i = 0; i < n; i++) {
			vector<double> vi(n);
			for (int j = 0; j < n; j++) {
				if (i == j) { vi[j] = 1.; }
				else { vi[j] = 0.; }
			}
			V.push_back(vi);
		}
	}
	// 检查S22中对角元为0的情况，若有，则返回其位置
	bool _check_diag(vector<vector<double>> S, const int r, const int l, int &pos) {
		const int n = S.size();
		bool check = false;
		for (int i = l; i < n - r - 1; i++) {
			if (S[i][i] == 0) check = true;
			pos = i;
		}
		return check;
	}
	// S22中有对角元S[pos][pos]为0，用Givens变换将第pos行全变为0，变换矩阵累乘到P上
	void _adapt(vector<vector<double>> &S, vector<vector<double>> &P,
		const int r, const int l, const int pos);
};

// make sure that m>=n
void SVD::decom(vector<vector<double>> A, vector<vector<double>> &P,
	vector<vector<double>> &S, vector<vector<double>> &Q)
{
	if (A.size() == 0 || A[0].size() == 0) { printf("err!\n"); return; }
	const int m = A.size(), n = A[0].size();
	if (m < n) { printf("本程序不处理m<n的情况，请将矩阵转置后重试\n"); return; }
	P.clear(); _init(P, m);
	S.clear();
	Q.clear(); _init(Q, n);

	// 二对角化 PAQ=S
	_hessenberg(P, A, Q); S = A;
	for (int i = n; i < m; i++) { S.pop_back(); }	// S是n*n矩阵
	
	// QR
	int iter = 0;	
	int r = 0, l = 0;				//每一步QR迭代的范围
	while (true) {					// 对S(H)完成schur分解
		iter++;
		_set_zero(S);				// 必要的置零操作
		_check_con(S, r, l);		// 检查是否收敛到对角阵，从上次继续
		///out::output(S,3); 
		///cout << "r & L " << r <<'\t'<< l << endl;
		if (r == n) { break; }
		int pos = 0;
		if (_check_diag(S, r, l, pos)) { _adapt(S, P, r, l, pos); }	// 检查是否有对角元为0
		else { _qr(P, S, Q, r, l); }			// 一次qr迭代		
	}
}

void SVD::_hessenberg(vector<vector<double>> &P,
	vector<vector<double>> &S, vector<vector<double>> &Q) {
	vector<double> v,x;
	double b;
	const int n = S[0].size(),m=S.size();
	for (int i = 0; i < n; i++) {
		x.clear();
		for (int k = i; k < m; k++) { x.push_back(S[k][i]); }

		funs::householder(x, v, b);
		funs::house_trans_l(S, v, b, i);
		funs::house_trans_l(P, v, b, i);

		if (i < n - 1) {
			x.clear();
			for (int k = i + 1; k < n; k++) { x.push_back(S[i][k]); }
			funs::householder(x, v, b);
			funs::house_trans_r(S, v, b, i+1);
			funs::house_trans_r(Q, v, b, i+1);
		}
	}
}

void SVD::_qr(vector<vector<double>>& P, vector<vector<double>>& S, vector<vector<double>>& Q, const int r, const int l)
{
	const int n = S.size()-r;// m=S[0].size()-r; // m=n
	if (n - l >= 3) {
		const double alpha = S[n - 1][n - 1] * S[n - 1][n - 1] + S[n - 2][n - 1] * S[n - 2][n - 1], beta = S[n - 2][n - 2] * S[n - 2][n - 1],
			delta = (S[n - 2][n - 2] * S[n - 2][n - 2] + S[n - 3][n - 2] * S[n - 3][n - 2] - alpha) / 2;
		const double mu = alpha - beta * beta / (delta + (delta >= 0 ? 1 : -1)*sqrt(delta*delta + beta * beta));
		double y, z,c,s;

		y = S[l][l] * S[l][l] - mu; z = S[l][l] * S[l][l + 1];
		for (int i = l; i < n-1; i++) {
			c = y / sqrt(y*y + z * z); s = -z / sqrt(y*y + z * z);
			funs::givens_trans_r(S, i, i + 1, c, s);
			funs::givens_trans_r(Q, i, i + 1, c, s);

			y = S[i][i]; z = S[i + 1][i];
			c = y / sqrt(y*y + z * z); s = z / sqrt(y*y + z * z);
			funs::givens_trans_l(S, i, i+1, c, s);
			funs::givens_trans_l(P, i, i+1, c, s);
			if (i < n - 2) {
				y = S[i][i + 1]; z = S[i][i + 2];
			}
		}
	}if (n - l == 2) {
		const double alpha = S[n-1][n-1] * S[n-1][n-1] + S[n-2][n-1] * S[n-2][n-1], beta = S[n-2][n-2] * S[n-2][n-1],
			delta = (S[n-2][n-2] * S[n-2][n-2] - alpha) / 2;
		const double mu = alpha - beta * beta / (delta + (delta >= 0 ? 1 : -1)*sqrt(delta*delta + beta * beta));
		double y, z, c, s;

		y = S[l][l] * S[l][l] - mu; z = S[l][l] * S[l][l + 1];
		c = y / sqrt(y*y + z * z); s = -z / sqrt(y*y + z * z);
		funs::givens_trans_r(S, l, l + 1, c, s);
		funs::givens_trans_r(Q, l, l + 1, c, s);

		y = S[l][l]; z = S[l + 1][l];
		c = y / sqrt(y*y + z * z); s = z / sqrt(y*y + z * z);
		funs::givens_trans_l(S, l, l + 1, c, s);
		funs::givens_trans_l(P, l, l + 1, c, s);
	}

}

inline void SVD::_set_zero(vector<vector<double>>& S)
{
	const int n = S.size();
	for (int i = 0; i < n - 1; i++) {
		if (fabs(S[i][i+1]) < u*(fabs(S[i][i]) + fabs(S[i + 1][i + 1]))) {
			S[i][i+1] = 0.;
		}
	}
	// 部分对角元置零
	double norm; funs::LInfnorm(S, norm);
	for (int i = 0; i < n; i++) {
		if (fabs(S[i][i]) < norm*u) { S[i][i] = 0.; }
	}
	/*
	// 次次对角元也进行置零操作 ，一般来说可做可不做，主要为了输出Schur阵时美观
	for (int i = 0; i < n - 2; i++) {
		for (int j = i + 2; j < n; j++) {
			if (fabs(S[j][i]) < u*fabs(S[i][i]) + u * fabs(S[j][j])) { S[j][i] = 0.; }
		}
	}*/
}

inline void SVD::_check_con(vector<vector<double>>& S, int & m, int & l)
{
	const int n = S.size();
	int i;
	for (i = n - m - 1; ; i--) {
		if (i <= 0) { m = n; break; }
		// 次对角元为0
		if (S[i-1][i] == 0) continue;
		/// H[i][i-1] != 0
		m = n - 1 - i; break;
	}// 右下角m*m矩阵是对角阵

	for (i = n - m - 1; i > 0; i--) {
		if (S[i-1][i] == 0) break;		// 再往下就不是不可约矩阵了
	}
	l = max(0, i); // 有 n-m-l != 1
	///cout << "L: " << l << endl;
}

void SVD::_adapt(vector<vector<double>>& S, vector<vector<double>>& P, const int r, const int l, const int pos)
{
	// S[pos][pos]=0 用Givens变换将整行变为0,即将S[pos][pos+1]变为1
	int n = S.size();
	double s, c;
	for (int i = pos + 1; i < n - r; i++) {
		s = -S[pos][i]/sqrt(S[pos][i]*S[pos][i]+S[i][i]*S[i][i]);
		c = S[i][i] / sqrt(S[pos][i] * S[pos][i] + S[i][i] * S[i][i]);
		funs::givens_trans_l(S, pos, i, c, s);	// 可提速，只对有变化的两列操作
		funs::givens_trans_l(P, pos, i, c, s);
	}
}
