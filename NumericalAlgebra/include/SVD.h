/* SVD �ֽ��㷨 */
/*
	�㷨���̣�
	���Խǻ�
	�������ж���QR����
		��r�ף����������Խ���ǰl�ף��������м�n-l-r������H22�ǵ����Ķ���
		QR��������H22�жԽ�ԪΪ0�������_adapt()����Ӧ��ȫ�����0
		  �������_qr()ִ�е���
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

	//�Ծ���A��SVD������A=P*S*Q������洢��������ά������
	void decom(vector<vector<double>> A, vector<vector<double>> &P, 
		vector<vector<double>> &S, vector<vector<double>> &Q);

private:
	double u;
	// ���Խǻ� P*S*Q = B
	void _hessenberg(vector<vector<double>> &P,
		vector<vector<double>> &S, vector<vector<double>> &Q);
	// ��S22��һ��qr�����������任�����۳˵�P��Q��
	void _qr(vector<vector<double>> &P, vector<vector<double>> &S,
		vector<vector<double>> &Q, const int r, const int l);
	// ����С�ķǶԽ�Ԫ����
	void _set_zero(vector<vector<double>> &S);
	// ���S���Խ������������
	void _check_con(vector<vector<double>> &S, int &m, int &l);
	// ��P��Q��ʼ��Ϊ��λ��
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
	// ���S22�жԽ�ԪΪ0����������У��򷵻���λ��
	bool _check_diag(vector<vector<double>> S, const int r, const int l, int &pos) {
		const int n = S.size();
		bool check = false;
		for (int i = l; i < n - r - 1; i++) {
			if (S[i][i] == 0) check = true;
			pos = i;
		}
		return check;
	}
	// S22���жԽ�ԪS[pos][pos]Ϊ0����Givens�任����pos��ȫ��Ϊ0���任�����۳˵�P��
	void _adapt(vector<vector<double>> &S, vector<vector<double>> &P,
		const int r, const int l, const int pos);
};

// make sure that m>=n
void SVD::decom(vector<vector<double>> A, vector<vector<double>> &P,
	vector<vector<double>> &S, vector<vector<double>> &Q)
{
	if (A.size() == 0 || A[0].size() == 0) { printf("err!\n"); return; }
	const int m = A.size(), n = A[0].size();
	if (m < n) { printf("�����򲻴���m<n��������뽫����ת�ú�����\n"); return; }
	P.clear(); _init(P, m);
	S.clear();
	Q.clear(); _init(Q, n);

	// ���Խǻ� PAQ=S
	_hessenberg(P, A, Q); S = A;
	for (int i = n; i < m; i++) { S.pop_back(); }	// S��n*n����
	
	// QR
	int iter = 0;	
	int r = 0, l = 0;				//ÿһ��QR�����ķ�Χ
	while (true) {					// ��S(H)���schur�ֽ�
		iter++;
		_set_zero(S);				// ��Ҫ���������
		_check_con(S, r, l);		// ����Ƿ��������Խ��󣬴��ϴμ���
		///out::output(S,3); 
		///cout << "r & L " << r <<'\t'<< l << endl;
		if (r == n) { break; }
		int pos = 0;
		if (_check_diag(S, r, l, pos)) { _adapt(S, P, r, l, pos); }	// ����Ƿ��жԽ�ԪΪ0
		else { _qr(P, S, Q, r, l); }			// һ��qr����		
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
	// ���ֶԽ�Ԫ����
	double norm; funs::LInfnorm(S, norm);
	for (int i = 0; i < n; i++) {
		if (fabs(S[i][i]) < norm*u) { S[i][i] = 0.; }
	}
	/*
	// �δζԽ�ԪҲ����������� ��һ����˵�����ɲ�������ҪΪ�����Schur��ʱ����
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
		// �ζԽ�ԪΪ0
		if (S[i-1][i] == 0) continue;
		/// H[i][i-1] != 0
		m = n - 1 - i; break;
	}// ���½�m*m�����ǶԽ���

	for (i = n - m - 1; i > 0; i--) {
		if (S[i-1][i] == 0) break;		// �����¾Ͳ��ǲ���Լ������
	}
	l = max(0, i); // �� n-m-l != 1
	///cout << "L: " << l << endl;
}

void SVD::_adapt(vector<vector<double>>& S, vector<vector<double>>& P, const int r, const int l, const int pos)
{
	// S[pos][pos]=0 ��Givens�任�����б�Ϊ0,����S[pos][pos+1]��Ϊ1
	int n = S.size();
	double s, c;
	for (int i = pos + 1; i < n - r; i++) {
		s = -S[pos][i]/sqrt(S[pos][i]*S[pos][i]+S[i][i]*S[i][i]);
		c = S[i][i] / sqrt(S[pos][i] * S[pos][i] + S[i][i] * S[i][i]);
		funs::givens_trans_l(S, pos, i, c, s);	// �����٣�ֻ���б仯�����в���
		funs::givens_trans_l(P, pos, i, c, s);
	}
}
