/* 展示作业成果 */
/* 
	exericse_1();	用过关Jacobi方法求实对称三对角阵的全部特征值和特征向量

	exercise_2();	二分法求实对称三对角阵特征值，再用反幂法求对应特征向量
*/

#include<vector>
#include<time.h>
#include"Solve.h"
#include"Symmetric.h"


void exercise_1() {
	cout << "--------------------- exercise 1 ----------------------" << endl;
	clock_t start_t, end_t;
	for (int n = 50; n <= 100; n+=10) {
		int iter;
		double time;
		vector<double> v;
		vector<vector<double>> A(n, vector<double>(n)),V;

		for (int i = 0; i < n - 1; i++) {
			A[i][i] = 4.;
			A[i][i + 1] = 1.;
			A[i + 1][i] = 1.;
		}
		A[n - 1][n - 1] = 4.;

		cout << "n = " << n << ":"<<endl;
		start_t = clock();
		Threshold_Jacobi J(1e-8);
		J.Jacobi(A, v, V, iter);
		end_t = clock();
		
		cout << "-------------------------------------------------------" << endl;
		cout << "this is Q:" << endl;
		out::output(V);
		cout << "-------------------------------------------------------" << endl;
		cout << "these are the eigenvalues:" << endl;
		sort(v.begin(), v.end());	// 从小到大排序
		out::output(v);
		cout << "iterate times: " << iter << endl;
		cout << "time consuming: "<<(end_t - start_t) / CLOCKS_PER_SEC << " s" << endl;
		cout << "-------------------------------------------------------" << endl;
	}
}

void exercise_2() {
	clock_t start_t, end_t;
	int n = 100;
	double eigenvalue, time;
	vector<double> v;
	vector<vector<double>> A(n, vector<double>(n));
	for (int i = 0; i < n-1; i++) {
		A[i][i] = 2.;
		A[i][i + 1] = -1.;
		A[i + 1][i] = -1.;
	}
	A[n - 1][n - 1] = 2.;

	divided D(1e-8);	//设置目标精度
	start_t = clock();
	D.eigensystem(A, 1, eigenvalue, v);
	end_t = clock();

	cout << "---------------------------------------------------" << endl;
	cout << "the minimum eigenvalue is: " << eigenvalue << endl;
	cout << "the eigen vector is: " << endl;
	out::output(v, 4);
	cout << "---------------------------------------------------" << endl;
	cout<< "time consuming: "<<(end_t - start_t) *1000/ CLOCKS_PER_SEC << " ms" << endl;
	cout << "-------------------------------------------------------" << endl;

	start_t = clock();
	D.eigensystem(A, 100, eigenvalue, v);
	end_t = clock();

	cout << "---------------------------------------------------" << endl;
	cout << "the maximum eigenvalue is: " << eigenvalue << endl;
	cout << "the eigen vector is: " << endl;
	out::output(v, 4);
	cout << "---------------------------------------------------" << endl;
	cout << "time consuming: " << (end_t - start_t)* 1000 / CLOCKS_PER_SEC  << " ms" << endl;
	cout << "-------------------------------------------------------" << endl;
}

int main() {
	//exercise_1();
	exercise_2();
}