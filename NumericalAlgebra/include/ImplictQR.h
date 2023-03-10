#pragma once
#include"Function.h"
#include<iostream>
#include<vector>

using namespace std;
// 用幂法求多项式模最大根，f是多项式系数（首项默认为1），max_iter是最大迭代次数，避免死循环
void Find_largest_root(vector<double> f, double &Largest_root, int &iter, double &times, const int max_iter=1000);

//
void Implict_qr(vector<vector<double>> A, vector<vector<double>> &schur);//vector<double> &lamb);

// A->H
void hessenberg_decom(vector<vector<double>> &A, vector<vector<double>> &V, vector<double> &b);

void set_zero(vector<vector<double>> &H, const double u);

void check_con(vector<vector<double>> &H, int &m, int &l);
// 
void two_step_qr(vector<vector<double>> &H, const int m, const int l);

void _get_sol_from_schur(vector<vector<double>> H, vector<vector<double>> &z);

void house_trans2(vector<vector<double>> &A, vector<double> v, const double beta, const int i_start);