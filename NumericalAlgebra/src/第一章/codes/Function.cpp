#include "Function.h"
#include <cmath>

using namespace std;

void forward_subs(vector<vector<double>>& L, vector<double>& b)
{
	// 前代法解方程 Lx = b，将结果存在变量b内，L是下三角阵，b是列向量 
	check(L.size(), b.size());
	for(int k=0; k<b.size(); k++){
		// 从前往后逐个求解 
		b[k] = b[k]/L[k][k];
		for(int j=k+1; j<b.size(); j++){
			b[j] = b[j] - b[k]*L[j][k];		// b[k] has been replaced by x[k]
		} 
	} 
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{
	// 前代法解方程 Lx = b，将结果存在变量b内，其中L是对角元为1的下三角阵，b是列向量，
	// 此时 L 可能不存储1，因此使用的数据范围在对角线以下（不含对角线） 
	check(L.size(), b.size());
	for(int k=0; k<b.size(); k++){
		for(int j=k+1; j<b.size(); j++){
			b[j] = b[j] - b[k]*L[j][k];
		} 
	} 
}

void back_subs(vector<vector<double>>& U, vector<double>& y)
{
	// 回代法解方程 Ux = y，结果存在向量b内，U是上三角阵，b是列向量
	check(U.size(), y.size());
	for(int k=y.size()-1; k>-1; k--){
		// 从后往前逐个求解 
		y[k] = y[k]/U[k][k];
		for(int j=k-1; j>-1; j--){
			y[j] = y[j] - y[k]*U[j][k];		// b[k] has been replaced by x[k]
		}
	}
}

void back_subs1(vector<vector<double>>& U, vector<double>& y) 
{
	// 回代法解方程 Ux = y，结果存在向量b内，U是对角元为1的上三角阵，b是列向量
	// 此时 U 可能不存储1，因此使用的数据范围在对角线以下（不含对角线） 
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
	// 高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
	// 参考教材上的算法
	for(int k=0; k<A.size(); k++){
		check(A[k][k]);
		// 计算L的第(k+1)列 
		for(int j=k+1; j<A.size(); j++){
			A[j][k] = A[j][k] / A[k][k];
		}
		// 让剩下部分适应前面的变换，自动计算U的第(k+1)行 
		for(int i=k+1; i<A.size(); i++){
			for(int j=k+1; j<A.size(); j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
			}
		}
	}
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{
	// 全主元高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
	// 参考教材上的算法，PAQ = LU  
	u.clear();
	v.clear();
	for(int k=0; k<A.size(); k++){
		
		// 选取主元 
		int p=k, q=k;		//记录主元位置 
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
		
		// 接下来和普通的高斯消元法一样 
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
	// 列主元高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
	// 参考教材上的算法，PA = LU 
	u.clear();
	for(int k=0; k<A.size(); k++){

		// 选取主元 
		int p=k;		//记录主元位置 
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

		// 接下来和普通的高斯消元法一样 
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
	// 行变换 P*b 
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
	// 列变换 Q*z
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
	// 对正定对称阵用cholesky算法分解，把结果储存在矩阵的下三角部分 
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
	
	// 补上上三角部分
	for(int i=0; i<A.size(); i++){
		for(int j=i+1; j<A.size(); j++){
			A[i][j] = A[j][i];
		}
	} 
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{
	// 对正定对称阵用改进的cholesky算法分解，把结果储存在矩阵的下三角部分 
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
	
	// 补上上三角部分
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
	//计算偏差 |Ax-b|，结果存在x中
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
	//检查矩阵行数与方程个数是否相符 
	if(sizeL != sizeb){
		printf("Error! Please check the number of equations.\n");
		exit(-1);
	}
}

void check(double a){
	//检查对角元是否为 0 
	if(fabs(a)<=0){//1e-16
		printf("Cannot solve the equation. Reason: the diagonal is too close to 0.\n");
		exit(1);
	}
}
