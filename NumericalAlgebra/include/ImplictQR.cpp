#include "ImplictQR.h"
#include<iostream>
#include<vector>
#include<Windows.h>
#include<stdlib.h>
#include<time.h>

//clock_t start_t, end_t;
using namespace std;

/*namespace*/
// 用幂法求多项式模最大根，f是多项式系数（次数从小到大，首项默认为1），max_iter是最大迭代次数，避免死循环
void Find_largest_root(vector<double> f, double &Largest_root, int &iter, double &times, const int max_iter) {
	// f <-> P
	const int n = f.size();
	double lamb, step;
	vector<double> u(n), y(n);//?
	for (int i = 0; i < n; i++) { u[i] = 1; } u[0] = 1;

	for (iter = 0; iter < max_iter; iter++) {
		// 先算出y[0]，后算yj，最后赋值y0（特殊的矩阵不调用已有函数）
		// 顺便记录其模最大值
		double tmp_y = 0; lamb = 0;
		for (int i = 0; i < n; i++) {// y0
			tmp_y -= f[n-1-i] * u[i];
		}
		for (int j = 1; j < n; j++) {
			y[j] = u[j - 1];
			if (fabs(lamb) < fabs(y[j])) { lamb = y[j]; }
		}
		y[0] = tmp_y; if (fabs(lamb) < abs(y[0])) { lamb = y[0]; }

		//求u，记录下u的变化幅度step
		step = 0.;
		for (int i = 0 ; i < n; i++) {
			if (fabs(u[i] - y[i] / lamb) > step) { step = abs(u[i] - y[i] / lamb); }
			u[i] = y[i] / lamb;
		}

		if (step < 1e-7) break;
		//test
		if (iter % 10 == 0){ cout << "iter:" << iter << "\tlambda:" << lamb << endl; }//
	}
	
	// 返回最大模根，注意lamb可能为负数
	Largest_root = lamb;
}


// 隐式QR算法求特征值，不改变A的值，但输出全体特征值lamb，加以修改可求特征向量Vec
void Implict_qr(vector<vector<double>> A, vector<vector<double>> &lamb) {
	vector<vector<double>> V; ///,Q
	vector<double> b;
	const int n = A.size();
	const double u = 1e-5;			//机器精度，后续改为传入式
	int m = 0, l = 0;				//每一步QR迭代的范围
	int iter = 0;

	///out::output(A);
	hessenberg_decom(A, V, b);		// A完成上Hessenberg分解，V和b是householder变换的向量集和比值集
	// hessenberg_decom(A, Q)
	///out::output(A);

	while (true) {			// 对A(H)完成schur分解
		iter++;
		set_zero(A, u);				// 必要的置零操作
		check_con(A, m, l);			// 检查是否收敛到schur阵，从上次继续
		//out::output(A,3); 
		if (m == n) { break; }
		///cout << "m & L " << m <<'\t'<< l << endl;
		two_step_qr(A, m, l);		// 一次双重步迭代
		// two_step_qr(A, m, l, Q); //同时记录Q	
	}

	out::output(A);
	printf("iterate times: %d\n", iter);
	_get_sol_from_schur(A, lamb);	// 从schur阵A中获取解
	///schur.clear();						// 先查看schur阵是否正确
	///for (int i = 0; i < n; i++) { schur.push_back(A[i]); }
}

// 上hessenberg分解，A->H，Householder变换的向量和系数分别存在V和b中
void hessenberg_decom(vector<vector<double>> &A, vector<vector<double>> &V, vector<double> &b) {
	V.clear();
	b.clear();
	const int n = A.size();
	for (int i = 0; i < n - 2; i++) {		// i< n-2
		vector<double> x(0), Vi(0);
		double bi;
		for (int j = i + 1; j < n; j++) { x.push_back(A[j][i]); }
		funs::householder(x, Vi, bi);		// 对第(i+1)列的第(i+2)~(n)分量应用Householder变换
		// x是n-1-i维的， Vi是n-2-i维的
		house_trans2(A, Vi, bi, i);

		V.push_back(Vi);
		b.push_back(bi);
	}

}

// 双重步位移的QR迭代，对其中的矩阵H22进行迭代，并将作用返回到H13和H23上
// 参考教材 P193
void two_step_qr(vector<vector<double>> &H, const int m, const int l) {// vector<vector<double>> &Q,
	const int n = H.size();
	///printf("hello?");

	double s = H[n-m-2][n-m-2] + H[n-m-1][n-m-1];
	double t = H[n-m-2][n-m-2] * H[n-m-1][n-m-1] - H[n-m-2][n-m-1] * H[n-m-1][n-m-2];
	double x = H[l][l] * H[l][l] + H[l][l+1] * H[l+1][l] - s * H[l][l] + t;
	double y = H[l+1][l] * (H[l][l] + H[l+1][l+1] - s);
	double z;	if (l+2 <= n-m-1) { z = H[l + 1][l] * H[l + 2][l + 1]; }
	///cout << "x,y,z " << x << '\t' << y << '\t'<< z << '\t' << endl;
	for (int k = l; k <= n - m - 3; k++) {
		vector<double> v, u = { x,y,z };
		double b;

		funs::householder(u, v, b);
		///out::output(v);
		// update(Q,m,l) (if nesseary)
		// 计算（q到n-1列）变换结果
		//int q = max(l+1, k);
		for (int j = 0; j < n; j++) {
			double beta;
			// beta = vT*bj; bj' = bj - beta * b * v
			// bj -> H[k:k+2][j]
			beta = H[k][j] + v[0] * H[k + 1][j] + v[1] * H[k + 2][j];
			H[k][j] -= beta * b;
			H[k + 1][j] -= beta * b*v[0];
			H[k + 2][j] -= beta * b*v[1];
		}
		// 计算（1到r行）变换结果
		int r = min(k + 3, n - m - 1);
		for (int i = 0; i <= r; i++) {
			double beta;
			beta = H[i][k] + H[i][k + 1] * v[0] + H[i][k + 2] * v[1];
			H[i][k] -= beta * b;
			H[i][k + 1] -= beta * b* v[0];
			H[i][k + 2] -= beta * b*v[1];
		}

		x = H[k + 1][k];
		y = H[k + 2][k];
		if (k < n - m - 3) { z = H[k + 3][k]; }
	}

	if (n - m - l == 2) {
		x = H[n - m - 2][n - m - 2];
		y = H[n - m - 1][n - m - 2];
	}	
	vector<double> v; double b;
	///cout << "x,y,z " << x << '\t' << y << endl;
	funs::householder({ x,y }, v, b);
	// update(Q,m,l) (if nesseary)
	
	// householder-每列的变换
	for (int j = max(0, n-m-3); j < n; j++) {
		double beta; 
		beta = H[n - m - 2][j] + v[0] * H[n - m - 1][j];
		H[n - m - 2][j] -= beta * b;
		H[n - m - 1][j] -= beta * b*v[0];
	}
	// householder-每行的变换
	for (int i = 0; i < n - m; i++) {
		double beta;
		beta = H[i][n - m - 2] + H[i][n - m - 1] * v[0];
		H[i][n - m - 2] -= beta * b;
		H[i][n - m - 1] -= beta * b*v[0];
	}
}

// 对上Hessenberg阵对角元的置零操作
void set_zero(vector<vector<double>> &H, const double u) {
	const int n = H.size();
	if (n != H[0].size()) { printf("err: size of H\n"); }

	for (int i = 0; i < n - 1; i++) {
		if (fabs(H[i+1][i]) < u*(fabs(H[i][i]) + fabs(H[i + 1][i + 1]))) {
			H[i+1][i] = 0.;
		}
	}
	// 次次对角元也进行置零操作 ，一般来说可做可不做，主要为了输出Schur阵时美观
	for (int i = 0; i < n - 2; i++) {
		for (int j = i + 2; j < n; j++) {
			if (fabs(H[j][i]) < u*fabs(H[i][i]) + u * fabs(H[j][j])) { H[j][i] = 0.; }
		}
	}
}

// 检查H到Schur阵的收敛进度
void check_con(vector<vector<double>> &H, int &m, int &l) {
	const int n = H.size();
	int i;
	for (i = n-m-1; ; i--) {
		if (i <= 0) { m = n; break; }
		// 次对角元为0
		if (H[i][i - 1] == 0) continue;
		/// H[i][i-1] != 0
		// 2*2 虚特征值子阵
		if ((H[i][i] - H[i-1][i-1])*(H[i][i] - H[i-1][i-1]) + H[i][i - 1] * H[i - 1][i] < 0) 
		{
			if (i < 2) { m = n; break; }
			if (H[i - 1][i - 2] == 0) { i--; continue; }
		}
		m = n - 1 - i; break; // H[i-1][i-2]!=0
	}// 右下角m*m矩阵是拟上三角阵

	for (i = n-m-1; i > 0; i--) {
		if (H[i][i - 1] == 0) break;		// 再往下就不是不可约矩阵了
	}
	l = max(0,i); // 有 n-m-l != 1
}

// 从Schur阵中提取特征值（可提取特征向量）
void _get_sol_from_schur(vector<vector<double>> H, vector<vector<double>> &z) {
	// x实部，y虚部
	const int n = H.size();
	z.clear();
	for (int i = 0; i < n-1; i++) {
		if (fabs(H[i + 1][i]) == 0) { z.push_back({ H[i][i], 0. }); }
		else {
			double x1 = (H[i][i] + H[i + 1][i + 1]) / 2;
			double delta = (H[i][i] - H[i + 1][i + 1])*(H[i][i] - H[i + 1][i + 1])/4 + H[i][i + 1] * H[i + 1][i];
			double y1 = sqrt(-delta);
			z.push_back({ x1, y1 });
			z.push_back({ x1, -y1 });
			i++;
		}
		if(i==n-2){ z.push_back({ H[i+1][i+1], 0. }); }
	}
	funs::check(z.size(), n);
}

void house_trans2(vector<vector<double>> &A, vector<double> v, const double beta, const int i_start) {
	// 计算A关于一轮householder变换后的相似矩阵 A2 = HT*A*H
	// 注意v是n-1-i_start维的，省略了第一个分量1
	const int n = A.size();
	const int m = v.size();
	if (m + i_start + 2 != n) { printf("Warning: house_trans2\n"); return; }//check

	// 计算分块 b = A[i_start+1:n-1, i_start:n-1] 变换后的结果，按列计算 HT*
	{
		for (int j = i_start; j < n; j++) {
			double k;
			// k = vT*bj
			k = A[i_start + 1][j];
			for (int i = 0; i < m; i++) { k += A[i_start + 2 + i][j] * v[i]; }
			// bj' = bj - beta k * v
			A[i_start + 1][j] -= beta * k;
			for (int i = 0; i < m; i++) { A[i_start + 2 + i][j] -= beta * k*v[i]; }
		}
	}
	// 计算分块 c = A[0:n-1, i_start+1:n-1] 变换后的结果，按行计算 *HT
	{
		for (int j = 0; j < n; j++) {
			double k;
			// k = bjT * v
			k = A[j][i_start + 1];
			for (int i = 0; i < m; i++) { k += A[j][i_start + 2 + i] * v[i]; }
			// bj' = bj - beta k * v
			A[j][i_start + 1] -= beta * k;
			for (int i = 0; i < m; i++) { A[j][i_start + 2 + i] -= beta * k*v[i]; }
		}
	}
}