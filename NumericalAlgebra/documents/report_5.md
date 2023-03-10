# 作业报告五
 - 陈泽高 PB20000302
## 代码结构：
* 本次作业含有四组代码，名称和功能如下
    * hw5.cpp：运行，直接展示作业结果
    * exercise.hpp：调用```Solve.hpp```中的求解器完成解方程任务，包括共轭梯度法（CG）、Jacobi迭代法、G-S迭代法以及SOR迭代法
    * Solve.cpp/hpp：以上一章作业为基础改编而成，实现了```Solver2```类和```SolverSparse```类，其中```Solver2```类是迭代法求解线性方程组的通用类，```SolverSparse```类是针对第一题的方程形式所编写的，```Solver```类是以前作业使用的类，本次作业没有用到它：
        1. ```Solver2```类的成员如下：
        ```
	    private:
		vector<vector<double>> A;	// 系数矩阵
		vector<double> b;		// 因变量
		vector<double> x;		// 解向量
		int iter;	    		// 最大迭代次数（设置为较大值，避免死循环）
        double tol;             // 迭代终止条件
        double w;               // 松弛因子
		
		iter_type decom;	    // 迭代方法
		
		double time_consuming = 0.0;    //单位：秒
		double loss = 0.0;		// 偏差：Ax-b的无穷范数
        double err;             // 当前迭代||x(k)-x(k+1)||
        ```
        方法：用Jacobi、G-S、SOR三种迭代格式求解线性方程组，同时记录下时间、偏差和迭代次数的方法。该部分功能调用Function.h的函数来实现。

        2. ```SolverSparse```类的成员如下：
        ```
        vector<vector<double>> g;	// 部分系数矩阵（即方程中的g和f）
        vector<vector<double>> f;	// 
        vector<vector<double>> x;	// 解矩阵u
        int iter = 1e5;				// 最大迭代次数
        double tol = 1e-7;			// 迭代终止条件
        double w = 1.2;				// 松弛因子

        double time_consuming = 0.0;// 耗时，单位：秒
        double err = 0.0;			// 当前迭代 ||x(k)-x(k+1)||Inf
        ```
        方法：与```Solver2```相似
    * Function.cpp/hpp：迭代法中等涉及矩阵运算的函数，打包到命名空间```funs```中，内容基本上是重新写的，调用的函数集中在后半部
	
	此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码；在```hw5.cpp```和```exercise.h```中可以调整部分参数，查看需要检查的解


## 一、二阶偏微分方程的差分（拉普拉斯方程）
* 问题描述：
偏微分方程
    $$ 
    \begin{cases}
    & -\Delta u + g(x,y) u = f(x,y), & (x,y) \in D \\
    & u(x,y) = \Gamma (x,y), & (x,y) \in \partial D\\
    & D = [0,1] \times [0,1]
    \end{cases}
    $$
同上题进行离散化后，同样得到大型稀疏方程组
    $$
    -u_{i-1,j} -u_{i,j-1} +(4+h^2 g(ih,jh))u_{i,j} -u_{i+1,j} -u_{i,j+1} = h^2 f(ih,jh)
    $$ $$1 \leq i,j \leq N-1$$
取 $g(x,y) = 1, f(x,y) = sin(x y), \Gamma (x,y) = x^2+y^2, n=20$ 进行计算
* 算法原理：
    1. 共轭梯度法：
    对于方程 $Ax=b$ ，选取了初值 $x_0$ 后，算出初始值 $r_0=b-A x_0$ ， $p_0 = r_0$ ， $r2_0=r_0^T r_0$，按照以下顺序迭代至达到收敛要求：
    $$
        q = Ap\\
        \alpha = r2 / p^T q\\
        r = r - \alpha q\\
        x = x + \alpha p\\
        tmp = r^T r\\
        \beta = tmp/r2\\
        r2 = tmp\\
        p = r + \beta p
    $$其中 $x$ 的步长为 $\alpha|p|_{\infty}$ ，可以在迭代过程中更新，本次作业所有采用共轭梯度法的代码均以此方式编写
    2. 本题的稀疏矩阵A具有特殊形式，可以针对该形式设计算法：
    对于 $x^T y$ 的计算，直接将二者作为向量拉直计算即可，详见文件```function.h```中的```dotSparse()```函数；
    对于 $q=A p$ 的计算，参考上次作业的算法可知：
    $$
    q(i,j) = (4 + h^2 g(i,j))p(i,j) - p(i-1,j) - p(i,j-1) - p(i+1,j) - p(i,j+1)
    $$其中 $A$ 在本题中由 $g(x,y)$ 确定，因此无需存储矩阵 $A$ ，直接通过上式计算即可，详见文件```function.cpp```文件中的函数```dotSparse()```
    3. 对于上述迭代的边界情况，可以对矩阵进行加边操作，上下左右各添加一行（列），其元素值均为 $0$ 。

* 运行结果及分析：
    1. 如图所示，共轭梯度法迭代次数 $54$ ，解得结果与SOR方法一致
    <center>
	<img src="screenshots/exe1_CG.png" style="zoom: 40%;" />

    <img src="screenshots/exe1_SOR.png" style="zoom: 40%;" />
    </center>
    
    2. 经过二分法筛选，得到最优迭代因子在 $1.3$ 左右，进一步细分得到的迭代次数仍然不低于 $62$
    <center>
	<img src="screenshots/exe1_SOR1.png" style="zoom: 40%;" /><img src="screenshots/exe1_SOR2.png" style="zoom: 40%;" /><img src="screenshots/exe1_SOR3.png" style="zoom: 40%;" /><img src="screenshots/exe1_SOR4.png" style="zoom: 40%;" />
    </center>

    不难看出，共轭梯度法在迭代次数方面占优，但是运算时间却略逊一筹。这主要是因为共轭梯度法一次迭代的计算量比SOR方法高出不少。

## 二、共轭梯度法求解方程组与其他方法对比
* 问题描述：
    1. 以Hilbert矩阵为系数矩阵，以所有自变量为 $\frac{1}{3}$ 为解的方程组
    2. 解方程$A x = b$，对比几种不同的迭代方法。其中
    $$ 
    A = \begin{pmatrix}
    10 &  1 &  2 &  3 &  4 \\
    1  &  9 & -1 &  2 & -3 \\
    2  & -1 &  7 &  3 & -5 \\
    3  &  2 &  3 & 12 & -1 \\
    4  & -3 & -5 & -1 & 15 \\
    \end{pmatrix}
    $$
    $$
    b = \begin{pmatrix}
    12&-27&14&-17&12 \\
    \end{pmatrix}^T
    $$

* 算法原理：
    参考前面提到的算法原理，略

* 运行结果与分析：
    1. 解Hilbert方程组，终止步长为 $10^{-7}$
    <center>
	<img src="screenshots/exe2.png" style="zoom: 50%;" />
    </center>

    可以看出，共轭梯度法能够快速准确地解出该方程组，求解时间以及精确度都远远超过之前编写的分解算法，也胜过Jacobi、G-S、SOR等迭代方法
    2. 解方程，对比不同迭代方法的时间和迭代次数
    <center>
	<img src="screenshots/exe3.png" style="zoom: 60%;" />
    </center>
    容易看出，共轭梯度法在运算时间和迭代次数方面都优于另外两个方法。由于该矩阵不是对角占优的，因此Jacobi迭代法和G-S迭代法都需要迭代更多次数才能达到好的收敛效果；对于五阶的方程组，由共轭梯度法的理论可知至多迭代五次，就可以达到极好的收敛效果，本次实验也印证了这一点。

