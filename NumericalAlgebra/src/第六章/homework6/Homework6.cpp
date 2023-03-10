/*
	本次作业主函数，用编写好的packages完成几个练习
	求方程的模最大根
		exercise1();
	求实矩阵全体特征值
		exercise2_1();//求高次多项式方程的所有根
		exercise2_2();//求矩阵的全部特征值，及其随其中某变元的变化
	用Mathematica计算相同的题目，参考结果在注释中，方便查阅
*/
#include "ImplictQR.h"

void exercise1() {
	vector<vector<double>> f;
	f.push_back({ 3., -5., 1. });			// x = -3
	f.push_back({ -1.,-3.,0. });			// x = 1.879385
	f.push_back({ -1000, 790, -99902, 79108.9, 9802.08, 10891.01, 208.01, 101 }); //x = -100

	for (int i = 1; i <= 3; i++) {
		cout << "------------------exercise 1."<< i <<"-------------------" << endl;
		double t, root;
		int iter;

		Find_largest_root(f[i-1], root, iter, t);			// 可添加参数：max_iter，表示最大迭代次数
		cout << "the Largest root is: " << root << endl;
		cout << "time consuming: " << ((t > 0.01) ? t : t * 1000)
			<< ((t > 0.01) ? "s" : "ms") << endl;
		cout << "iterate times: " << iter << endl;
	}
}

void exercise2_1() {
	cout << "------------------exercise 2.1-------------------" << endl;
	int N = 41; //矩阵大小 
	int iter;

	vector<vector<double>> V, D(N, vector<double>(N));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			D[i][j] = 0.;
		}
	}
	for (int i = 0; i < N - 1; i++) {
		D[i + 1][i] = 1.;
	}
	D[0][N - 4] = -1.; D[0][N - 1] = -1.;
	Implict_qr(D, V, iter);
	cout << "iterate times: " << iter << endl;
	cout << "方程的解如下（每行一个解，" 
		<< "每个解由实部（左）及虚部（右）组成，"
		<< "若无虚部则为实根）:" << endl;
	out::output(V);
	/*
		reference solutions:
		-0.952484,
		-0.96814 - 0.120867 I,
		-0.96814 + 0.120867 I,
		-0.956339 - 0.273776 I,
		-0.956339 + 0.273776 I,
		-0.910511 - 0.425528 I, 
		-0.910511 + 0.425528 I,
		-0.836863 - 0.567826 I,
		-0.836863 + 0.567826 I,
		-0.739101 - 0.695904 I,
		-0.739101 + 0.695904 I,
		-0.620673 - 0.805889 I,
		-0.620673 + 0.805889 I,
		-0.48528 - 0.894538 I,
		-0.48528 + 0.894538 I, 
		-0.336984 - 0.959228 I,
		-0.336984 + 0.959228 I,
		-0.180206 - 0.997962 I,
		-0.180206 + 0.997962 I,
		-0.0197286 - 1.00935 I,
		-0.0197286 + 1.00935 I,
		0.139165 - 0.992477 I,
		0.139165 + 0.992477 I,
		0.289812 - 0.946424 I,
		0.289812 + 0.946424 I,
		0.417152 - 0.871067 I,
		0.417152 + 0.871067 I, 
		0.507569 - 0.810573 I,
		0.507569 + 0.810573 I,
		0.63234 - 0.753401 I,
		0.63234 + 0.753401 I,
		0.75372 - 0.65538 I,
		0.75372 + 0.65538 I, 
		0.855158 - 0.532634 I,
		0.855158 + 0.532634 I, 
		0.933664 - 0.392546 I,
		0.933664 + 0.392546 I,
		0.987184 - 0.240354 I,
		0.987184 + 0.240354 I, 
		1.0143 - 0.080923 I,
		1.0143 + 0.080923 I
	*/
}

void exercise2_2() {
	cout << "------------------exercise 2.2-------------------" << endl;
	vector<vector<double>> A = { {9.1, 3.0, 2.6, 4.0}, {4.2, 5.3, 4.7, 1.6},
		{3.2, 1.7, 9.4, 0.}, {6.1, 4.9, 3.5, 6.2} }, Eigen;
	// A[2][3]待定
	int iter;

	for (double x = 0.9; x < 1.2; x += 0.1) {
		A[2][3] = x;
		Implict_qr(A, Eigen, iter);
		cout << "When x = " << x << ":" << endl;
		cout << "iterate times: " << iter << endl;
		cout << "The eigenvalues are: " << endl;
		int n = Eigen.size();
		for (int i = 0; i < n; i++) {
			if (Eigen[i].size() > 1) {
				cout << Eigen[i][0] << " + " << Eigen[i][1] << "I"
					<< '\t' << "模长: "
					<< sqrt(Eigen[i][0] * Eigen[i][0] + Eigen[i][1] * Eigen[i][1])
					<< endl;
			}
			else {
				cout << Eigen[i][0] << '\t' << '\t'
					<< '\t' << "模长: "
					<< fabs(Eigen[i][0])
					<< endl;
			}
		}
	}
}

int main() {

	exercise1();

	exercise2_1();

	exercise2_2();
}
