#include "exercise.h"
//#include

int main() {
/**/
	{//exercise 1
		const int n = 100;
		const double a = 0.5;
		const double tol = 1e-6;

		double epsilon = 1;
		//exercise_1(n, a, epsilon, tol);
		//system("pause");
		epsilon = 0.1;
		exercise_1(n, a, epsilon, tol);
		system("pause");
		epsilon = 0.01;
		exercise_1(n, a, epsilon, tol);
		system("pause");
		epsilon = 0.0001;
		exercise_1(n, a, epsilon, tol);
		system("pause");
	}
/**/
	{//exercise 2
		int N;
		N = 20;
		exercise_2(N);
		system("pause");
		N = 40;
		exercise_2(N);
		system("pause");
		N = 60;
		exercise_2(N);
	}
	//cout << fabs(0) << endl;
}
