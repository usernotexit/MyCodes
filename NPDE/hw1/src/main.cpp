#include<stdio.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<Eigen/Core>

using namespace std;

// Warning: you should make sure that Func returns a double variable
template<typename Func>
void getValue(Func f, size_t m, Eigen::VectorXd& fx)
{
    fx.resize(m);
    double deltaX = 2 * EIGEN_PI / (double)m;
    for (int i = 0; i < m; i++)
    {
        fx[i] = f(i * deltaX);
    }
}

void printValue(const Eigen::VectorXd& fx)
{
    //cout << fixed << setprecision(3);
    //cout.setf(ios::scientific);
    auto n = fx.size();
    for (int i = 0; i < n; i++)
    {
        cout << fx[i] << '\t';
    }
    cout << endl;
}

void fourier(int N, size_t m)
{
    // functions
    auto xx = [&](double const x)->double {return x; };
    auto v = [&](double const x)->double {return (EIGEN_PI-x)/2; };
    auto vN = [&](double const x)->double {
        double f = 0;
        for (int j = 1; j <= N; j++)
            f += sin(j * x) / j;
        return f;
    };
    auto vtN = [&](double const x)->double {
        double f = 0;
        for (int j = 1; j <= N; j++)
            f += sin(j * x) / j * sin(j*EIGEN_PI/N) / (j*EIGEN_PI/N);
        return f;
    };

    // values
    Eigen::VectorXd xj, Vj, VNj, VtNj, V_VNj, V_VtNj;
    getValue(xx, m, xj);
    getValue(v, m, Vj);
    getValue(vN, m, VNj);
    getValue(vtN, m, VtNj);
    V_VNj = Vj - VNj;
    V_VtNj = Vj - VtNj;

    // print out v(x)-vN(x)
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "m = " << m << endl;
    cout << "v(x) - vN(x): " << endl;
    printValue(V_VNj);
    cout << "v(x) - vtN(x): " << endl;
    printValue(V_VtNj);

    // print v, vN, vtN on txt, then on excel
    ofstream file;
    file.open("Fourier.txt", std::ios::app|ios::out);
    if(!file)
    {
        std::cerr << "Unable to open file!" << endl;
        return;
    }

    file << "x" << '\t' << "v(x)" << '\t' << "vN(x)" << '\t' << "v(x)-vN(x)"
        << '\t' << "vtN(x)" << '\t' << "v(x)-vtN(x)" << endl;
    for (int j = 0; j < m; j++)
    {
        file << xj[j] << '\t' << Vj[j] << '\t' << VNj[j] << '\t' << V_VNj[j]
            << '\t' << VtNj[j] << '\t' << V_VtNj[j] << endl;
    }
    file.close();
}

int main()
{
    fourier(10, 20);
    fourier(100, 160);

    return 1;
}