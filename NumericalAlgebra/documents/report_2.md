# 作业报告二
## 陈泽高 PB20000302
### 代码结构：
* 本次作业含有四组代码，名称和功能如下
    * Homework2.cpp：运行终端，直接展示作业结果
    * Solve.cpp/hpp：以上一章作业为基础改编而成，实现了```Solver```类，包含如下内容：
        1. 成员函数：
        ```
	    private:
		vector<vector<double>> A;	// 系数矩阵
		vector<vector<double>> A_;	// 分解后矩阵
		vector<double> b;		// 因变量
		vector<double> x;		// 解向量
		vector<int> u;			// 行变换矩阵
		vector<int> v;			// 列变换矩阵
		
		decom_type decom = not_yet;	// 分解标识
		LU_type lu = none;
		
		double time_consuming = 0.0;    //单位：秒
		double loss = 0.0;		// 偏差：Ax-b的无穷范数
		double err = 0.0;	        // 估计相对误差上界
        ```
        2. 用gauss消元法和Cholesky分解法求解线性方程组，同时记录下时间和偏差的方法。该部分功能调用Function.h的函数来实现。此外，还将**矩阵分解**与**解上（下）三角线性方程组**的步骤分开，便于单独对系数矩阵进行分析
        3. 分析矩阵$A$的条件数的方法和分析方程$Ax=b$解的相对误差上界的方法，该部分功能调用Err.h的函数实现
    * Function.cpp/hpp：Gauss消元法、Cholesky分解法等涉及矩阵运算的函数、还有检查函数check()等，打包到命名空间```funs```中，内容与第一章作业基本一致
    * Err.h: 针对第二章作业编写的扩展，编写了计算矩阵和向量的无穷范数的算法（直接求法和优化法）
	
	此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码。部分在本次作业中不需要使用的函数没有写完，对作业没有影响。本次实现的类可以移植到其他程序中，也可以根据新的需求进行扩展。

### 一、估算条件数
* 问题描述：估算5到30阶Hilbert矩阵的条件数
* 算法原理：条件数的定义
$$
\kappa_{\infty}(A)=||A||_{\infty}*||A^{-1}||_{\infty}
$$
    1. 其中$||A||_{\infty}$容易计算，选择绝对值行和最大者即可
    2. 对于$||A^{-1}||_{\infty}$，可以用优化法计算，具体算法参考教材2.5.1，其中求$w$和$z$的步骤由原本的矩阵乘法转化为解线性方程组，在系数矩阵A已经完成LU分解的条件下，算法的时间复杂度为$O(n^2)$，核心代码如下（伪代码）：
    ```
    // A是系数矩阵，AT是A的转置，n是矩阵行列数，x、w、v、w、z等都是临时使用的变量，此处Solve返回解向量
	for(i=0; i<n; i++){ x=1/n; }
	bool is_min = false;
	for (i = 0; i < n; i++) {   // 转置
		for (j = 0; j < n; j++) {
			AT[i][j] = A[j][i];
		}
	}

	while(!is_min){
	    // w = (A^-1)T x; v = sign(w)
	    w = Solve(AT, x);
	    for (i = 0; i < n; i++) {
	    	v[i] = sign(w[i]);
	    }
    
	    // z = (A^-1)*v
	    z = Solve(A, v);

	    // ||z|| vs zT*x
	    normInf_z = 0.0, zx=0.0;
	    LInfnorm(z, normInf_z, k);  //详见源代码,normInf_z即为z的无穷范数，k指示z绝对值最大元素的位置
	    for(i=0; i<n; i++){ zx += x[i]*z[i]; }
		
	    if (normInf_z <= fabs(zx)) { 
            norm = normInf_z;
            is_min = true; 
        }	// areadly minimized
		
        else{
            x = 0;      // x的所有元素变成0
            x[k] = 1;
	    }
	}
    ```
    由上述算法得到的是$norm = ||A^{-1}||_{\infty}$

* 运行结果与分析：
    <center>

	<img src="ex1.jpg" style="zoom: 70%;" />
	
	***exercise 1***
	</center>

    上图显示了5到30阶Hilbert矩阵条件数的估算值，当n较小时与参考值基本一致，随着n的增大，与参考值开始出现偏差。此外，一开始笔者采用了列主元法对A进行LU分解，最终得到的条件数和普通消元法的结果并不相同，该问题可能与Hilbert矩阵的病态性质相关。（在Homework.cpp中调整相应的注释即可查阅）

### 二、估算解的精度
* 问题描述：对于给定的一个系数矩阵$A$（见教材），随机选取解$x\in \mathbb{R}^n$，算出因变量$b=Ax$，然后解方程$Ax=b$，估计解的精度下界，比较机器解$x_1$和真实解$x$之间的相对误差
* 算法描述：课堂上推导过解的偏差上界
    $$
    \frac{||x-x_1||_{\infty}}{||x||_{\infty}} = 
    \kappa_{\infty}(A) \frac{||r||_{\infty}}{||b||_{\infty}}
    $$
    其中$r = Ax_1 - b$，有前面求解条件数的算法，这一部分不难解决，详见源代码（Err.h中的函数```Solver::error()```）
    需要注意的是，该题给出的系数矩阵的条件数较小，因此需要选取较大的$x$，且最好不为整数，才能让求出的误差不会太小。
* 运行结果与分析：
    <center>

	<img src="ex2.jpg" style="zoom: 70%;" />
	
	***exercise 2***
	</center>
    
    其中$loss = ||r||$，$err$为估计误差（即上式右项），$true err$是实际误差（即上式左项）
    可以看出，该矩阵的条件数较小，得到的解误差小，不过随着矩阵阶数的增加，估计误差和实际误差不可避免地增长。
    总体来看，实际误差会比估计误差小一个数量级左右。