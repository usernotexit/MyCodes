#include<iomanip>
#include "exercise.h"

using namespace std;

int main() {
	/*
	{//exe1
		// f=sin(x*y), phi=x^2+y^2, g=1
		const int n = 20;
		const double tol = 1e-7, h = 1./n;
		vector<vector<double>> g(n + 1, vector<double>(n + 1)),x;
		vector<vector<double>> f(n + 1, vector<double>(n + 1));
		vector<vector<double>> phi(n + 1, vector<double>(n + 1));
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				g[i][j] = h * h;
				f[i][j] = h * h * sin(i*h*j*h);
			}
			phi[0][i] = i * h * i * h;//edge
			phi[n][i] = 1 + phi[0][i];
			phi[i][0] = phi[0][i];
			phi[i][n] = phi[n][i];
		}

		exe1_CG(n, g, f, phi, tol);
		for (int i = 0; i < 10; i++) {
			double w = 1.725 + 0.001*i;
			exe1_SOR(n, g, f, phi, tol, w);//1.7~1.8 //1.7~1.75
		}
	}*/
	/*
	{//exe2
		int n = 20;
		exe2_Hilbert(n);
		n = 40;
		exe2_Hilbert(n);
		n = 60;
		exe2_Hilbert(n);
		n = 80;
		exe2_Hilbert(n);
	}*/
	/**/
	{//exe3
		vector<vector<double>>A =
		{ {10, 1,2,3,4},
		{1,9,-1,2,-3},
		{2,-1,7,3,-5},
		{3,2,3,12,-1},
		{4,-3,-5,-1,15} };
		vector<double> b = { 12,-27,14,-17,12 };

		Solver2 s(A, b, 1e-6, 1000);
		exe3_solve(s, CG);
		exe3_solve(s, Jacobi);
		exe3_solve(s, GS);
	}
}
