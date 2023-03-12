//#include <cstdio>
#include <chrono>
#include <iostream>
#include <vector>

#include "Eigen/Eigen"
//#include <Eigen/QR>

bool Assert(bool state);

// todo 5: change the class in to a template class
//#define Tem double
template<typename Tem>
class Matrix {
private:
	int rows, cols;
	Tem *data;

public:
    // default constructor
    // https://en.cppreference.com/w/cpp/language/constructor
	Matrix() :rows(1), cols(1), data(new Tem[1]) { data[0] = 0; }	// Matrix() = delete;

    // constructor with initilizer list
	Matrix(int r, int c)
		: rows(r), cols(c), data(new Tem[r*c]) {}		// can we add default 0 ?
	Matrix(const Matrix& Mat)							// 引用和常值Mat都可以使用该重载
		: rows(Mat.rows), cols(Mat.cols), data(new Tem[Mat.rows * Mat.cols])
	{
		// std::cout << "Created!" << std::endl;
		for (int n = 0; n < rows*cols; n++)
			data[n] = Mat[n];
	}

    // desctructor
    // https://en.cppreference.com/w/cpp/language/destructor
	~Matrix() { if(data!=nullptr)delete data; }

    int nrow() const {return rows;}
    int ncol() const {return cols;}

    // operator overloding
    Tem& operator()(int r, int c) {return data[r*cols + c]; }		// 可能改变矩阵的“引用”
	Tem& operator()(int r, int c) const { return data[r*cols + c]; }// 不改变矩阵值的“引用”
	Tem& operator[](int n) { return data[n]; /* todo 3: particular entry of the matrix*/ }
    Tem& operator[](int n) const {return data[n]; /* todo 3: particular entry of the matrix*/}

    Matrix col(int col) 
	{/* todo 4: particular column of the matrix*/
		Matrix Mat_Col(rows, 1);			// 会不会在该函数结束时被析构？不会，因为返回的是一个值，不是引用
		for (int i = 0; i < rows; i++)
			Mat_Col[i] = data[i*cols+col];
		return Mat_Col;
	}
    Matrix row(int row) 
	{/* todo 4: particular row of the matrix*/
		Matrix Mat_Row = Matrix(1, cols);
		for (int i = 0; i < cols; i++)
			Mat_Row(0, i) = (*this)(row,i);// data[row*cols + i];
		return Mat_Row;
	}
    Matrix submat(int rowL, int rowU, int colL, int colU) const 
	{/* todo 4: return a sub-matrix specified by the input parameters*/
		Assert((rowL<=rowU) && (colL<=colU));
		Matrix Mat_Sub(rowU-rowL+1, colU-colL+1);
		for (int i = 0; i < rowU - rowL + 1; i++)
		{
			for (int j = 0; j < colU - colL + 1; j++)
				Mat_Sub(i, j) = (*this)(i + rowL, j + colL);//data[(i + rowL)*cols + (j + colL)];
		}
		return Mat_Sub;
	}

    // constant alias
	Matrix& operator= (const Matrix& rhs)
	{
		// std::cout << "'=' called\n ~~~~~~~~~~~" << std::endl;
		if(data!=nullptr) delete[] data;
		Tem* newData = new Tem[rhs.cols*rhs.rows];
		rows = rhs.rows; cols = rhs.cols;
		for (int i = 0; i < cols*rows; i++)
			newData[i] = rhs.data[i];

		data = newData;
		return *this;
	}

	/* 这里采取引用为返回值的好处在于方便后续 +=、-= 等运算符的重载 
	弊端则是对于 (A*B).print() 之类的命令会出错，因为返回的是局部变量的引用 */
    Matrix operator+ (const Matrix& rhs) 
	{
		Assert(cols == rhs.ncol() && rows == rhs.nrow());
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] + rhs[n];
		return res;
	}
    Matrix operator- (const Matrix& rhs) 
	{
		Assert(cols == rhs.ncol() && rows == rhs.nrow());
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] - rhs[n];
		return res;
	}
	
	// 矩阵乘法
    Matrix operator* (const Matrix& rhs) 
	{
		Assert(cols == rhs.nrow());
		Matrix res(rows, rhs.ncol());
		for(int i=0; i<rows; i++)
			for (int j = 0; j < rhs.ncol(); j++)
			{
				res(i, j) = 0;
				for (int k = 0; k < cols; k++)
					res(i, j) += (*this)(i, k) * rhs(k, j);
			}
		// res.print(); //test
		return res;
	}

 	// 对应位置元素相除，请不要把带有0元素的矩阵放在除号右侧
	Matrix operator/ (const Matrix& rhs) 
	{
		Assert(rows == rhs.nrow() && cols == rhs.ncol());
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] + rhs[n];
		return res;
	}

    Matrix operator+= (const Matrix& rhs) 
	{
		Matrix res(*this + rhs);
		*this = res;
		return *this;
	}

    Matrix operator-= (const Matrix& rhs) 
	{
		Matrix res(*this - rhs);
		*this = res;
		return *this;
	}

    Matrix operator*= (const Matrix& rhs) 
	{
		Matrix res(*this * rhs);
		*this = res;
		return *this;
	}
	Matrix operator/= (const Matrix& rhs)
	{
		Matrix res(*this / rhs);
		*this = res;
		return *this;
	}

	Matrix operator+ (Tem v)
	{
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] + v;
		return res;
	}
    Matrix operator- (Tem v) 
	{
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] - v;
		return res;
	}
    Matrix operator* (Tem v) 
	{
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] * v;
		return res;
	}
    Matrix operator/ (Tem v) 
	{
		Assert(fabs(v) < 1e-10);	// 要求v不要太接近0
		Matrix res(rows, cols);
		for (int n = 0; n < cols*rows; n++)
			res[n] = (*this)[n] / v;
		return res;
	}

	Matrix operator+= (Tem v)
	{
		Matrix res(*this + v);
		*this = res;
		return *this;
	}
	Matrix operator-= (Tem v)
	{
		Matrix res(*this - v);
		*this = res;
		return *this;
	}
    Matrix operator*= (Tem v) 
	{
		Matrix res(*this * v);
		*this = res;
		return *this;
	}
    Matrix operator/= (Tem v) 
	{
		Matrix res(*this / v);
		*this = res;
		return *this;
	}

    void print () const {
        printf("this matrix has size (%d x %d)\n", rows, cols);
        printf("the entries are:\n");
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
				std::cout << (*this)(i, j) << '\t';
			std::cout << std::endl;
		}
    }
	
};

bool Assert(bool state)
{
	if (state)
		return true;
	else
		std::cout << "Err!" << std::endl;
	//exit(0);
}


// BONUS: write a sparse matrix class
template<typename R>
class SparseMatrix {

};

//template<typename type>
//void test()
int main() 
{
	/* test 1.1: 准确性测试 */
    Matrix<int> A(3,4), B(4,3), R(4,3);			// 指定行列的建构函数
	srand(time(0));

	for (int i = 0; i < A.nrow(); i++)			// 方法 nrow(), ncol()
		for (int j = 0; j < A.ncol(); j++)
			A(i, j) = i + j;					// 重载操作符 '()'
	for (int i = 0; i < B.nrow(); i++)
		for (int j = 0; j < B.ncol(); j++)
		{
			B(i, j) = i - j;
			R(i, j) = rand() % 5;
		}
	
	Matrix<int> C = A*B;						// 构造函数 Matrix(const Matrix& Mat)
	Matrix<int> D; D = B*A;						// 默认建构函数，重载运算符 '='

	(B+R).print();								// 重载运算符 '*' 要求返回Matrix值
    
	A.print();									// 函数 print() ，输出矩阵值
    B.print();
    C.print();
	D.print();
	
	/* test 1.2 准确性测试 */
	// https://en.cppreference.com/w/cpp/chrono
	auto start = std::chrono::steady_clock::now();	
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> runTime = end - start;
	std::cout << "Runtime: " << runTime.count() << std::endl;
    
	// todo 8: use Eigen and compare
	//auto D = Eigen::Array22d(1, 34, 3, 4);
	//std::cout << D;
	Eigen::Matrix<double, 3, 4>();

    return 0;
}
