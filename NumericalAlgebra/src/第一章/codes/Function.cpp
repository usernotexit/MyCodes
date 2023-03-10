#include "Function.h"
#include <cmath>

using namespace std;

void forward_subs(vector<vector<double>>& L, vector<double>& b)
{
	// ǰ�����ⷽ�� Lx = b����������ڱ���b�ڣ�L����������b�������� 
	check(L.size(), b.size());
	for(int k=0; k<b.size(); k++){
		// ��ǰ���������� 
		b[k] = b[k]/L[k][k];
		for(int j=k+1; j<b.size(); j++){
			b[j] = b[j] - b[k]*L[j][k];		// b[k] has been replaced by x[k]
		} 
	} 
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{
	// ǰ�����ⷽ�� Lx = b����������ڱ���b�ڣ�����L�ǶԽ�ԪΪ1����������b����������
	// ��ʱ L ���ܲ��洢1�����ʹ�õ����ݷ�Χ�ڶԽ������£������Խ��ߣ� 
	check(L.size(), b.size());
	for(int k=0; k<b.size(); k++){
		for(int j=k+1; j<b.size(); j++){
			b[j] = b[j] - b[k]*L[j][k];
		} 
	} 
}

void back_subs(vector<vector<double>>& U, vector<double>& y)
{
	// �ش����ⷽ�� Ux = y�������������b�ڣ�U����������b��������
	check(U.size(), y.size());
	for(int k=y.size()-1; k>-1; k--){
		// �Ӻ���ǰ������ 
		y[k] = y[k]/U[k][k];
		for(int j=k-1; j>-1; j--){
			y[j] = y[j] - y[k]*U[j][k];		// b[k] has been replaced by x[k]
		}
	}
}

void back_subs1(vector<vector<double>>& U, vector<double>& y) 
{
	// �ش����ⷽ�� Ux = y�������������b�ڣ�U�ǶԽ�ԪΪ1����������b��������
	// ��ʱ U ���ܲ��洢1�����ʹ�õ����ݷ�Χ�ڶԽ������£������Խ��ߣ� 
	check(U.size(), y.size());
	for(int k=y.size()-1; k>-1; k--){
		y[k] = y[k]/U[k][k];
		for(int j=k-1; j>-1; j--){
			y[j] = y[j] - y[k]*U[j][k];		// b[k] has been replaced by x[k]
		}
	}
}

void gauss_elim(vector<vector<double>>& A)
{
	// ��˹��ȥ������A����LU�ֽ⣬���ֱ�洢��A���ϰ벿�֣����Խ��ߣ����°벿�֣������Խ��ߣ� 
	// �ο��̲��ϵ��㷨
	for(int k=0; k<A.size(); k++){
		check(A[k][k]);
		// ����L�ĵ�(k+1)�� 
		for(int j=k+1; j<A.size(); j++){
			A[j][k] = A[j][k] / A[k][k];
		}
		// ��ʣ�²�����Ӧǰ��ı任���Զ�����U�ĵ�(k+1)�� 
		for(int i=k+1; i<A.size(); i++){
			for(int j=k+1; j<A.size(); j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
			}
		}
	}
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{
	// ȫ��Ԫ��˹��ȥ������A����LU�ֽ⣬���ֱ�洢��A���ϰ벿�֣����Խ��ߣ����°벿�֣������Խ��ߣ� 
	// �ο��̲��ϵ��㷨��PAQ = LU  
	u.clear();
	v.clear();
	for(int k=0; k<A.size(); k++){
		
		// ѡȡ��Ԫ 
		int p=k, q=k;		//��¼��Ԫλ�� 
		for(int i=k; i<A.size(); i++){
			for(int j=k; j<A.size(); j++){
				if(fabs(A[i][j]) > fabs(A[p][q])) {
					p = i;
					q = j;
				}
			}
		}
		for(int i=0; i<A.size(); i++){	// exchange row k with row p 
			double swamp; 
			swamp = A[k][i]; 
			A[k][i] = A[p][i]; 
			A[p][i] = swamp;
		}
		for(int j=0; j<A.size(); j++){	//exchange column k with column q
			double swamp; 
			swamp = A[j][k]; 
			A[j][k] = A[j][q]; 
			A[j][q] = swamp;
		}
		u.push_back(p);
		v.push_back(q);
		
		// ����������ͨ�ĸ�˹��Ԫ��һ�� 
		check(A[k][k]);
		for(int j=k+1; j<A.size(); j++){
			A[j][k] = A[j][k] / A[k][k];
		}
		for(int i=k+1; i<A.size(); i++){
			for(int j=k+1; j<A.size(); j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
			}
		}
	}
}

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
{
	// ����Ԫ��˹��ȥ������A����LU�ֽ⣬���ֱ�洢��A���ϰ벿�֣����Խ��ߣ����°벿�֣������Խ��ߣ� 
	// �ο��̲��ϵ��㷨��PA = LU 
	u.clear();
	for(int k=0; k<A.size(); k++){

		// ѡȡ��Ԫ 
		int p=k;		//��¼��Ԫλ�� 
		for(int i=k; i<A.size(); i++){
			if(fabs(A[i][k]) > fabs(A[p][k])){
				p = i;
			}
		}
		for(int i=0; i<A.size(); i++){	// exchange row k with row p 
			double swamp;
			swamp = A[k][i];
			A[k][i] = A[p][i];
			A[p][i] = swamp;
		}

		u.push_back(p);

		// ����������ͨ�ĸ�˹��Ԫ��һ�� 
		check(A[k][k]);
		for(int j=k+1; j<A.size(); j++){
			A[j][k] = A[j][k]/A[k][k];
		}
		for(int i=k+1; i<A.size(); i++){
			for (int j=k+1;  j<A.size(); j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
			}
		}
	}
}

void vector_pb(vector<int>& u, vector<double>& b)
{
	// �б任 P*b 
	// u represents matrix P
	
	// check(u, b.size());
	for(int k=0; k<b.size(); k++){
		double swap;
		swap = b[k];
		b[k] = b[u[k]];
		b[u[k]] = swap;
	}
}

void vector_qb(vector<int>& v, vector<double>& b)
{
	// �б任 Q*z
	// v represents matrix Q
	
	// check(v, b.size());
	for(int k=b.size()-1; k>=0; k--){	// the order is just reversed
		double swap;
		swap = b[k];
		b[k] = b[v[k]];
		b[v[k]] = swap;		
	}
}

void cholesky_decomp(vector<vector<double>>& A)
{
	// �������Գ�����cholesky�㷨�ֽ⣬�ѽ�������ھ���������ǲ��� 
	for(int k=0; k<A.size(); k++){
		A[k][k] = sqrt(A[k][k]);
		for(int j=k+1; j<A.size(); j++){
			A[j][k] = A[j][k]/A[k][k];
		}
		for(int j=k+1; j<A.size(); j++){
			for(int i=j; i<A.size(); i++){
				A[i][j] = A[i][j] - A[i][k]*A[j][k];
			}
		}
	}
	
	// ���������ǲ���
	for(int i=0; i<A.size(); i++){
		for(int j=i+1; j<A.size(); j++){
			A[i][j] = A[j][i];
		}
	} 
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{
	// �������Գ����øĽ���cholesky�㷨�ֽ⣬�ѽ�������ھ���������ǲ��� 
	for(int j=0; j<A.size(); j++){
		std::vector<double> v(A.size());
		for(int i=0; i<j; i++){
			v[i] = A[j][i]*A[i][i];
			A[j][j] = A[j][j] - A[j][i]*v[i];
		}
		for(int i=j+1; i<A.size(); i++){
			for(int k=0; k<j; k++){
				A[i][j] = A[i][j] - A[i][k]*v[k];
			}
			A[i][j] = A[i][j]/A[j][j];
		}
	}
	
	// ���������ǲ���
	for(int i=0; i<A.size(); i++){
		for(int j=i+1; j<A.size(); j++){
			A[i][j] = A[j][i];
		}
	} 
}

void matrix_DLT(vector<vector<double>>& A)
{
	for(int i=0; i<A.size(); i++){
		for(int j=i+1; j<A.size(); j++){
			A[i][j] = A[i][j]*A[i][i];
		}
	}
}

void vector_Ax(vector<vector<double>> A, vector<double> x, vector<double> b, double &detla) {
	//����ƫ�� |Ax-b|���������x��
	detla = 0.0;
	check(A.size(), x.size());
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			b[i] -= A[i][j] * x[j];
		}
		detla += b[i] * b[i];
	}
	detla = sqrt(detla);
}


void check(int sizeL, int sizeb){
	//�����������뷽�̸����Ƿ���� 
	if(sizeL != sizeb){
		printf("Error! Please check the number of equations.\n");
		exit(-1);
	}
}

void check(double a){
	//���Խ�Ԫ�Ƿ�Ϊ 0 
	if(fabs(a)<=0){//1e-16
		printf("Cannot solve the equation. Reason: the diagonal is too close to 0.\n");
		exit(1);
	}
}
