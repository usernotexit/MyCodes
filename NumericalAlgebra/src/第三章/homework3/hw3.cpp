#include "QR.h"
#include "Err.h"

using namespace std;
void solve_all(vector<vector<double>> A, vector<double> b, int N) {
	vector<double> x; double loss, time;
	// vector<Solver> S;
	Solver S1(A, b), S2(A, b), S3(A, b), S4(A, b);
	cout << "QR method..." << endl;
	S1.QR();
	S1.solve();
	S1.showshow_result(x, time, loss);
	for (int i = 0; i < N; i++) { cout << x[i] << '\t'; }
	cout << endl;
	cout << "loss: " << loss << endl;
	cout << "time comsuming: " << time << "s" << endl;

	cout << "gauss eliminate..." << endl;
	S2.gauss_elim();
	S2.solve();
	S2.showshow_result(x, time, loss);
	for (int i = 0; i < N; i++) { cout << x[i] << '\t'; }
	cout << endl;
	cout << "loss: " << loss << endl;
	cout << "time comsuming: " << time << "s" << endl;

	cout << "gauss eliminate column pivoting..." << endl;
	S3.gauss_elim_col_pivoting();
	S3.solve();
	S3.showshow_result(x, time, loss);
	for (int i = 0; i < N; i++) { cout << x[i] << '\t'; }
	cout << endl;
	cout << "loss: " << loss << endl;
	cout << "time comsuming: " << time << "s" << endl;

	cout << "gauss eliminate full pivoting..." << endl;
	S4.gauss_elim_full_pivoting();
	S4.solve();
	S4.showshow_result(x, time, loss);
	for (int i = 0; i < N; i++) { cout << x[i] << '\t'; }
	cout << endl;
	cout << "loss: " << loss << endl;
	cout << "time comsuming: " << time << "s" << endl;
}
void solve_cholesky(vector<vector<double>> A, vector<double> b, int N) {
	vector<double> x; double loss, time;
	Solver S1(A, b), S2(A, b);
	cout << "cholesky..." << endl;
	S1.cholesky_decomp();
	S1.solve();
	S1.showshow_result(x, time, loss);
	for (int i = 0; i < N; i++) { cout << x[i] << '\t'; }
	cout << endl;
	cout << "loss: " << loss << endl;
	cout << "time comsuming: " << time << "s" << endl;

	cout << "modified cholesky..." << endl;
	S2.modified_cholesky_decomp();
	S2.solve();
	S2.showshow_result(x, time, loss);
	for (int i = 0; i < N; i++) { cout << x[i] << '\t'; }
	cout << endl;
	cout << "loss: " << loss << endl;
	cout << "time comsuming: " << time << "s" << endl;
}

int main() {
	/*	*/
	{// exercise 1.1
		int N = 84; //¾ØÕó´óÐ¡

		// initialize
		std::vector<vector<double>> A(N, vector<double>(N));
		std::vector<double> b(N);
		for (int i = 0; i < N - 1; i++)
		{
			A[i][i] = 6;
			A[i + 1][i] = 8;
			A[i][i + 1] = 1;
			b[i] = 15;
		}
		A[N - 1][N - 1] = 6;
		b[0] = 7;
		b[N - 1] = 14;
		solve_all(A, b, N);
	system("pause");
	}

	{// exercise 1.2
		//initialize
		int N = 100;
		std::vector<std::vector<double>> A(N, vector<double>(N));
		std::vector<double> b(N);
		for (int i = 0; i < N - 1; i++) {
			A[i][i] = 100;
			A[i][i + 1] = 1;
			A[i + 1][i] = 1;
		}
		A[N - 1][N - 1] = 100;
		for (int i = 0; i < N; i++) {
			b[i] = (double)rand() / (double)RAND_MAX * 100;
		}
		solve_all(A, b, N);
		solve_cholesky(A, b, N);
		system("pause");
	}

	
	{// exe 1.3
		//initialize
		int N = 40;
		std::vector<vector<double>> A(N, std::vector<double>(N));
		std::vector<double> b(N);
		for (int i = 0; i < N; i++) {
			b[i] = 0;
			for (int j = 0; j < N; j++) {
				A[i][j] = 1.0 / (i + j + 1);
				b[i] += 1.0 / (i + j + 1);
			}
		}
		solve_all(A, b, N);
		solve_cholesky(A, b, N);
		system("pause");
	}
	{// exe 2
		vector<vector<double>> A,B;
		vector<double> t = { -1,-0.75,-0.5,0,0.25,0.5,0.75 };
		vector<double> y = { 1,0.8125,0.75,1,1.3125,1.75,2.3125 };
		vector<double> abc;
		double min;
		for (int i = 0; i < t.size(); i++) {
			vector<double> a = { t[i]*t[i], t[i], 1 };
			A.push_back(a);
		}
		Solver S(A, y);
		cout << "Solving LS problem..." << endl;
		S.solve_LS(abc, min);
		cout << "result(a,b,c): " << abc[0] << '\t' << abc[1] << '\t' << abc[2] << endl;	
		cout << "loss: " << min << endl;
		system("pause");
	}

	{// exe 3
		vector<vector<double>> A =
		{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
		{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
		{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
		{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
		{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
		{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
		{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
		{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
		{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
		{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
		{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
		{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
		{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
		{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
		{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
		{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
		{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
		{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
		{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
		{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
		{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
		{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
		{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
		{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
		{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
		{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
		{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
		{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
		vector<double> b =
		{ 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
		28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
		30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
		37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 };
		double min;
		vector<double> x;
		Solver S(A, b);
		cout << "Solving LS problem..." << endl;
		S.solve_LS(x, min);
		for (int i = 0; i < x.size(); i++) { cout << x[i] << '\t'; }
		cout << endl;
		cout << "loss: " << min << endl;

		system("pause");
	}
}

