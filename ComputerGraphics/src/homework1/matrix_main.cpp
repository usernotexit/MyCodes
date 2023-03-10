#include <cstdio>

constexpr auto T = double;
// todo 5: change the class in to a template class
class Matrix {
public:
    // default constructor
    // https://en.cppreference.com/w/cpp/language/constructor
	delete Matrix() { }

    // constructor with initilizer list
	Matrix(int a, int b)
		: rows(a), cols(b), data(new T[a*b](0)) {}		// can we add default 0 ?

    // desctructor
    // https://en.cppreference.com/w/cpp/language/destructor
	~Matrix() { delete[] data; }

    int nrow() const {return rows;}
    int ncol() const {return cols;}

    // operator overloding
    double operator()(int r, int c) {return &data[r + rows*c]; }
    double operator[](int) {return &data[n]; /* todo 3: particular entry of the matrix*/}

    Matrix col(int) {return Matrix(); /* todo 4: particular column of the matrix*/}
    Matrix row(int) {return Matrix(); /* todo 4: particular row of the matrix*/}
    Matrix submat(int, int, int, int) const {return Matrix(); /* todo 4: return a sub-matrix specified by the input parameters*/}

    // constant alias
    Matrix& operator= (const Matrix& rhs) {/*.todo 3..*/return *this;}

    Matrix operator+ (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}
    Matrix operator- (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}
    Matrix operator* (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}
    Matrix operator/ (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}

    Matrix operator+= (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}
    Matrix operator-= (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}
    Matrix operator*= (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}
    Matrix operator/= (const Matrix& rhs) {/*.todo 3..*/ return Matrix();}

    Matrix operator+ (double v) {/*.todo 3..*/return Matrix();}
    Matrix operator- (double v) {/*.todo 3..*/return Matrix();}
    Matrix operator* (double v) {/*.todo 3..*/return Matrix();}
    Matrix operator/ (double v) {/*.todo 3..*/return Matrix();}

    Matrix operator+= (double v) {/*.todo 3..*/return Matrix();}
    Matrix operator-= (double v) {/*.todo 3..*/return Matrix();}
    Matrix operator*= (double v) {/*.todo 3..*/return Matrix();}
    Matrix operator/= (double v) {/*.todo 3..*/return Matrix();}

    void print () const {
        printf("this matrix has size (%d x %d)\n", rows, cols);
        printf("the entries are:\n");
        /* todo 4: print all the entries of the matrix */
    }

private:
int rows, cols;
double *data;
};

// BONUS: write a sparse matrix class
template<typename R>
class SparseMatrix {

}

int main() {
    Matrix A, B;

    // todo 6: fill A anb B with random numbers
    //for(int i=0; i<A.nrow(); i++)
    //    for(int j=0; j<A.ncol(); j++)
    //        A(i, j) = rand();

    // todo 7: benchmark with runtime, using std::chrono
    // https://en.cppreference.com/w/cpp/chrono


    Matrix C = A*B;
    A.print();
    B.print();
    C.print();

    // todo 8: use Eigen and compare

    return 0;
}
