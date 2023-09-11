#pragma once
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

enum decom_type { not_yet, Gauss, Gauss_ful, Gauss_col, Cholesky, modified_Cholesky, QR };
enum LU_type { none, L, U };
enum iter_type { Jacobi, GS, SOR, CG };

namespace funs {
	// ��������
	void LInfnorm(vector<vector<double>> &A, double &norm);// ����A�������
	void LInfnorm(vector<vector<double>> &A, vector<vector<double>> &B, double &norm);// ����A-B�������
	void LInfnorm(const vector<double> b, double &norm, int &j);//����b�������
	void LInfnorm(const vector<double> b, double &norm);//����b�������
	void LInfnorm(const vector<double> a, const vector<double> b, double &norm);// a-b �������
	void loss_Inf(vector<vector<double>> A, vector<double> x, vector<double> b, double &detla);//����ƫ�� Ax-b
	void loss_2(vector<vector<double>> A, vector<double> x, vector<double> b, double &loss);// ����ƫ�� Ax-b
	void dot(vector<double> x, vector<double> y, double &result);//����x��y���ڻ�
	void times(vector<vector<double>> A, vector<double> x, vector<double> &y);//�������˷�A*x
	double dot(vector<double> x, vector<double> y);//���������ڻ���ֱ�ӷ��ؽ��
	void minus(vector<double> x, vector<double> y, vector<double> &z);// ����x-y=z
	void dot(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C); // ����˷�C=A B
	void minus(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C); // �������C=A-B
	void transpose(vector<vector<double>> &A); //����ת��
	// householder
	void householder(vector<double> x, vector<double> &v, double &b);//����Householder�ֽ⣬��һ��
	void house_trans_l(vector<vector<double>> &A, vector<double> v, double b, int i_start);//�������A��householder�任����ˣ�
	void house_trans_r(vector<vector<double>> &A, vector<double> v, double b, int j_start);//�������A��householder�任���ҳˣ�
	// givens
	void givens_trans_r(vector<vector<double>> &A, const int p, const int q, const double c, const double s);// �ҳ�Givens
	void givens_trans_l(vector<vector<double>> &A, const int p, const int q, const double c, const double s);// ���Givens

	// ������ �ⷽ��
	void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����
	void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����
	void back_subs(vector<vector<double>>& U, vector<double>& y);//�ش���
	void back_subs1(vector<vector<double>>& U, vector<double>& y);//�Խ�ԪΪ1�Ļش���
	void vector_pb(vector<int>&u, vector<double>&b);//��������P*b����ѡ��
	void vector_qb(vector<int>& v, vector<double>& b);//��������Q*b����ѡ��
	void matrix_DLT(vector<vector<double>>& A);//�������D*L^T����ѡ��
	// �ֽⷨ
	void gauss_elim(vector<vector<double>>& A);//Gauss��ȥ��
	void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//ȫ��ԪGauss��ȥ��
	void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��
	void cholesky_decomp(vector<vector<double>>& A);//�Գ��������׼Cholesky�ֽ�
	void modified_cholesky_decomp(vector<vector<double>>& A);//�Ľ���ƽ������
	// QR�ֽ�
	void QR_decom(vector<vector<double>> &A, vector<double> &v);	//QR�ֽ����
	void house_trans(vector<vector<double>> &A, vector<double> v, const double b, const int i, const int j);//����Ծ���A��householder�任
	void house_trans(vector<double> &b, vector<double> v, const double beta, const int j_start);//�����b��househoer�任��ע������龳��Ҫ��
	// ������
	void updateJacobi(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y);// ����Jacobi������һ�θ���
	void updateGS(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y, double w);// ����G-S��SOR������һ�θ���
	void updateJacobiSparse(const vector<vector<double>> g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y);//����Jacobi��������һ�θ���
	void updateGSSparse(const vector<vector<double>> g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y);//����G-S��������һ�θ���
	void updateSORSparse(const vector<vector<double>> g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y, const double w);//����SOR��������һ�θ���
	void updateCGSparse(const vector<vector<double>>g, const vector<vector<double>> f, vector<vector<double>> x, vector<vector<double>> &y );//���㹲���ݶȷ���һ�θ���
	double dot(vector<vector<double>> x, vector<vector<double>> y);//�����ڻ�x*y��x��y��ֱ
	void timesSparse(vector<vector<double>> g, vector<vector<double>> p, vector<vector<double>> &q);//�������˷�A*p��A��g����
	double dotSparse(vector<vector<double>> x, vector<vector<double>> y);//�����ڻ�x*y��x��y��ֱ // specially for Sparse

	void check(int sizeL, int sizeb);	//�����������뷽�̸����Ƿ���� 
	void check(double a);				//���Խ�Ԫ�Ƿ�Ϊ 0 

}


namespace out {
	// outputϵ�к���ר������������������չʾ
	void output(vector<double> a, int pre=6);

	void output(vector<vector<double>> A, int pre=6);
}
