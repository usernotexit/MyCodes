# 作业报告三
## 陈泽高 PB20000302
### 代码结构：
* 本次作业含有四组代码，名称和功能如下
    * Homework3.cpp：运行终端，直接展示作业结果
    * Solve.cpp/hpp：以上一章作业为基础改编而成，实现了```Solver```类，包含如下内容：
        1. 本次作业使用的成员如下：
        ```
	    private:
		vector<vector<double>> A;	// 系数矩阵
		vector<vector<double>> A_;	// 分解后矩阵
		vector<double> b;		// 因变量
		vector<double> x;		// 解向量
		vector<int> u;			// 行变换矩阵
		vector<int> v;			// 列变换矩阵
        vector<int> beta;       // Householder分解的β值
		
		decom_type decom = not_yet;	// 分解标识
		LU_type lu = none;
		
		double time_consuming = 0.0;    //单位：秒
		double loss = 0.0;		// 偏差：Ax-b的无穷范数
        ```
        2. 用gauss消元法、Cholesky分解法和QR分解求解线性方程组，同时记录下时间和偏差的方法。该部分功能调用Function.h和QR.h的函数来实现。此外，还将**矩阵分解**与**解上（下）三角线性方程组**的步骤分开，便于单独对系数矩阵进行分析
    * Function.cpp/hpp：Gauss消元法、Cholesky分解法等涉及矩阵运算的函数，打包到命名空间```funs```中，内容在第一章作业的基础上，增加了实现Householder变换和QR分解的功能函数```householder()```、```house_trans()```、```QR_decom()```，调整了部分函数，以适配系数矩阵不是方阵时求解LS问题
    * QR.hpp：针对第三章作业编写```Solver```类的拓展，包括系数矩阵的QR分解和求解LS问题的方法```Solver::QR()```，```Solver::solve_LS()```
	
	此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码
### 一、利用QR分解求解线性方程组
* 问题描述：利用QR分解方法分解系数矩阵，再求解线性方程组
* 算法原理：
    1. 利用Householder变换的性质：
        * 不改变向量的$2$范数（即正交变换）
        * 对任意向量$\pmb{x}$，可以构造Householder变换$H(\pmb{v}) = \pmb{I}-2\beta \pmb{v} \pmb{v}^T$，使得
        $$
        H(\pmb{v}) \pmb{x} = \alpha \pmb{e}_1
        $$
        详见教材算法或者源代码，核心代码如下：
    ```
    void householder(vector<double> x, vector<double> &v, double &b) {
	    // 输入向量x，计算householder变换Q使得Qx只有第一个分量不为0，其中Q用v和b表示，详见教材算法
	    // Q = I - b v*vT
	    // v的第一分量固定为1，不保留
	    LInfnorm(x, normInf);//x的L无穷范数
	    v0_tmp = x[0] / normInf;
	    x[0] = v0_tmp;
	    for (int i = 1; i < n; i++) { x[i] = x[i] / normInf; v.push_back(x[i]); }
	    halfnorm = 0;
	    for (int i = 1; i < x.size(); i++) { halfnorm += x[i] * x[i]; }

	    if (halfnorm <= 1e-16 ) { b = 0; }
	    else 
	    {
	    	a = sqrt(x[0] * x[0] + halfnorm);
	    	if (x[0] <= 0) { v0_tmp = x[0] - a; }
	    	else { v0_tmp = -halfnorm / (x[0] + a); }
	    	b = 2 * v0_tmp * v0_tmp / (halfnorm + v0_tmp * v0_tmp);
	    	for (unsigned int i = 1; i < n; i++) { v[i - 1] = v[i - 1] / v0_tmp; }	// n-1 elements in vector v
	    }
	}
    ```
    2. 对系数矩阵逐列进行Householder变换，将其正交变换为上三角阵，即QR分解,核心代码如下：
    ```
    void QR_decom(vector<vector<double>> &A, vector<double> &d) {
	    // 对行数为m，列数为n的矩阵A进行QR分解，结果储存在A的上三角部分和下半部分，每一次变换的beta存在数组d中
	    // 注意使用前面两条函数前。要先进行适应性调整（例如增加起始位置(i,j)以便在原矩阵A上截取子矩阵）
	    int m = A.size(), n = A[0].size();
	    for (int j = 0; j < n; j++) {
	        double b;
	        vector<double> x, v;
	        for (int i = j; i < m; i++) { x.push_back(A[i][j]); }

	        householder(x, v, b);
	        house_trans(A, v, b, j, j); // 对A进行Householder变换

	        d.push_back(b);
	        for (int i = j + 1; i < m; i++) {
	            A[i][j] = v[i - j - 1];
	        }
        }
    }
    ```
    3. 分解之后，可以求解线性方程组，核心代码如下：
    ```
    void Solver::solve(){
	    funs::check(A_.size(), b.size());
	    x.clear();
	    for (unsigned int i = 0; i < A_[0].size(); i++) { x.push_back(b[i]); }

	    // QR
	    if (decom == 6) {
	    	start_t = clock();
		    int m = A_.size(), n = A_[0].size();
		    if (m < n) { cout << "We can't solve this now, try new data" << endl; return; }

		    // caculate QT*b
		    for (int i = 0; i < A_[0].size(); i++) {
		        vector<double> v;
		        v.push_back(1.0);
		        for (int j = i + 1; j < int(A_.size()); j++) {
                    v.push_back(A_[j][i]);
		        }
		        funs::house_trans(x, v, beta[i], i);
		    }
        }
    }
    ```
* 运行结果与分析：
    1. 用gauss、Cholesky和QR方法分别求解方程组，比较各种方法优劣

	<img src="screenshots/p1_1.jpg" style="zoom: 40%;" />
    <center>

	***exercise 1.1***
    </center>
    <center>
	<img src="screenshots/s1_1_qr.jpg" style="zoom: 50%;" />
    </center>   
    <center>
	<img src="screenshots/s1_1_gauss.jpg" style="zoom: 50%;" />
	</center>   
    <center>
	<img src="screenshots/s1_1_gausscol.jpg" style="zoom: 50%;" />
	</center>   
    <center>
	<img src="screenshots/s1_1_gaussfull.jpg" style="zoom: 50%;" />
	</center>    

    <div STYLE="page-break-after: always;"></div>
	<img src="screenshots/p1_2.jpg" style="zoom: 40%;" />
    <center>

	***exercise 1.2***
    </center>
	<center>
	<img src="screenshots/s1_2_qr.jpg" style="zoom: 50%;" />
	</center>   
    <center>
	<img src="screenshots/s1_2_gauss.jpg" style="zoom: 50%;" />
	</center>   
    <center>
	<img src="screenshots/s1_2_gausscol.jpg" style="zoom: 50%;" />
	</center>   
    <center>
	<img src="screenshots/s1_2_gaussfull.jpg" style="zoom: 50%;" />
	</center>   
    <center>
    <img src="screenshots/s1_2_cholesky.jpg" style="zoom: 50%;" />
	</center>   

    <div STYLE="page-break-after: always;"></div>
    <center>

	***exercise 1.3***
    </center>
	以四十阶Hilbert矩阵为系数，解全为1的方程
	<center>
	<img src="screenshots/s1_3_qr.jpg" style="zoom: 50%;" />
	</center>   
    <center>
	<img src="screenshots/s1_3_gauss.jpg" style="zoom: 50%;" />
	</center>   
    <center>
    <img src="screenshots/s1_3_cholesky.jpg" style="zoom: 50%;" />
	</center> 

    可以看出，在上述几例当中，QR分解方法的速度略快于高斯分解，慢于cholesky分解；解的精度略差于列主元高斯分解

### 二、利用QR分解求解LS问题
* 问题描述：利用QR分解方法分解系数矩阵，再求解最小二乘问题
* 算法原理：

    4. 求解LS问题：QR分解已经完成，可以直接得到LS问题的最优解。核心函数如下：
    ```
    void Solver::solve_LS() {
        if (decom == 0) { this->QR(); }     // 对系数矩阵A分解，得到正交阵Q和R
        if (decom != 6) { cout << "please do QR decomposition" << endl; return; }

        int m = A_.size(), n = A[0].size();

        for (int i = 0; i < b.size(); i++) { x.push_back(b[i]); }

        // caculate QT*b
        for (int i = 0; i < A_[0].size(); i++) {
            vector<double> v;
            v.push_back(1.0);           // 补上前面不保留的1
            for (int j = i + 1; j < int(A_.size()); j++) {
                v.push_back(A_[j][i]);
            }
            funs::house_trans(x, v, beta[i], i);
        }

        // get (c1, c2)T = QT*b
        for (int i = m; i > n; i--) { x.pop_back(); }

        // solve Rx = c1;
        funs::back_subs(A_, x);
    }
    ```

* 运行结果与分析：
    1. 求二次多项式$y = a t^2 + b t + c$，使得残向量的$2$范数最小的意义下拟合表中数据
    <center>
    <img src="screenshots/p2.jpg" style="zoom: 60%;" />
	</center> 
    <center>
    <img src="screenshots/s2.jpg" style="zoom: 60%;" />
	</center> 
    较好地完成拟合

    2. 求房地产估计的线性模型
    <center>
    <img src="screenshots/p3.jpg" style="zoom: 50%;" />
	</center>     
    <center>
    <img src="screenshots/s3.jpg" style="zoom: 60%;" />
	</center> 
    产生了一定偏差，该模型的准确性有待进一步验证
