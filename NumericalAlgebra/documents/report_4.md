# 作业报告四
## 陈泽高 PB20000302
### 代码结构：
* 本次作业含有四组代码，名称和功能如下
    * hw4.cpp：运行终端，直接展示作业结果
    * exercise.hpp：利用```Solve.hpp```中的求解器完成解方程任务，分别采用Jacobi迭代法、G-S迭代法以及SOR迭代法
    * Solve.cpp/hpp：以上一章作业为基础改编而成，实现了```Solver2```类和```SolverSparse```类，其中```Solver2```类是迭代法求解线性方程组的通用类，```SolverSparse```类是针对第二题的方程形式所编写的，```Solver```类是以前作业使用的类，本次作业没有用到它：
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
	
	此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码；由于篇幅原因笔者没有把方程的具体解放入报告，在```exercise.cpp```中可以调整输出，查看需要检查的解

### 一、求差分方程的解
* 问题描述：
两点边值问题
    $$ 
    \left\{
    \begin{aligned}
    & \epsilon \frac{\mathrm{d}^2 y}{\mathrm{d} x^2} + \frac{\mathrm{d}y}{\mathrm{d}x}= a, 0<a<1 \\
    & y(0) = 0 , y(1) = 1\\
    \end{aligned}
    \right.
    $$  
    将区间$[0,1]$ $n$等分，对方程离散化后，原方程近似为大型线性方程组$A y = b$，其中
    $$ 
    A = \begin{pmatrix}
    -(2 \epsilon + h) & \epsilon + h \\
    \epsilon & -(2 \epsilon + h) & \epsilon + h \\
    & \epsilon & -(2 \epsilon + h) & \ddots \\
    &  & \ddots & \ddots & \epsilon + h\\
    &  &  & \epsilon & -(2 \epsilon + h) \\
    \end{pmatrix}
    $$
    $$
    y = \begin{pmatrix}
    y_{1} & y_{2} & \cdots & y_{n-1} \\
    \end{pmatrix}^T
    $$
    $$
    b = \begin{pmatrix}
    a h^2 & a h^2 & \cdots & a h^2 & a h^2 - (\epsilon + h) \\
    \end{pmatrix}^T
    $$
    目标：求解线性方程组并比较与微分方程精确解的偏差
* 算法原理：
参考教材给出的迭代算法，每次迭代的代码如下（代码已经非常直观，不多作赘述）：
    1. Jacobi迭代法
    ```
    void updateJacobi(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y) {
        // y = ((L+U)*x + b) * D^-1
        // x是更新前的向量，y是更新后的向量
        y.clear();

        for (int i = 0; i < A.size(); i++) {
            double tmp = 0;
            for (int j = 0; j < A[0].size(); j++) {
                if (j == i)continue;
                else tmp -= A[i][j] * x[j];
            }
            tmp += b[i];
            y.push_back(tmp / A[i][i]);
        }
    }
    ```
    2. G-S迭代法和SOR迭代法（二者写在一起，G-S法中$w=1$）
        ```
        void updateGS(const vector<vector<double>> A, const vector<double> b, vector<double> x, vector<double> &y, double w) {
            // D * y = (1-w)x + w * (L*y + U*x + b)
            y.clear();
            for (int i = 0; i < x.size(); i++) { y.push_back(x[i]); }

            for (int i = 0; i < A.size(); i++) {
                double tmp = 0;
                for (int j = 0; j < A[0].size(); j++) {
                    if (j == i) continue;
                    else tmp -= A[i][j] * y[j];
                }
                tmp += b[i];
                y[i] = (1 - w) * y[i] + w * tmp / A[i][i];
            }
        }
        ```
* 运行结果与分析：
    1. $\epsilon = 1$ 时部分输出截图（最后一张的时间过长是由于中途暂停导致的）：
    <center>
	<img src="screenshots/1_1_J.jpg" style="zoom: 50%;" />
    <img src="screenshots/1_1_G.jpg" style="zoom: 50%;" />

    <img src="screenshots/1_1_SOR1.jpg" style="zoom: 50%;" /><img src="screenshots/1_1_SOR2.jpg" style="zoom: 50%;" /><img src="screenshots/1_1_SOR3.jpg" style="zoom: 50%;" />
    </center>

    2. 三种算法的收敛性符合预期，且得到的解十分相近，以下为 $\epsilon = 0.1$ 时不同算法给出的解
    <center>
	<img src="screenshots/1_2_J.jpg" style="zoom: 40%;" /><img src="screenshots/1_2_G.jpg" style="zoom: 40%;" />
    <img src="screenshots/1_2_SOR1.jpg" style="zoom: 40%;" /><img src="screenshots/1_2_SOR2.jpg" style="zoom: 40%;" /><img src="screenshots/1_2_SOR3.jpg" style="zoom: 40%;" />
    </center>

    3. $\epsilon$ 不同时，利用不同迭代法的迭代次数、求解时间、偏差如下：

        | 迭代次数 | Jacobi | G-S | SOR w=1.1 | w=1.2 | w=1.3 |
        | :---: | :----: | :----: | :----: | :----: | :----: |
        | $\epsilon = 1$ | 15403 | 6576 | 5547 | 4660 | 3882 |
        | $\epsilon = 0.1$ | 6050 | 2980 | 2496 | 2082 | 1724 |
        | $\epsilon = 0.01$ | 569 | 332 | 280 | 235 | 195 |
        | $\epsilon = 0.0001$ | 117 | 108 | 141 | 不收敛 | 不收敛 |

        | 耗时（单位：s） | Jacobi | G-S | SOR w=1.1 | w=1.2 | w=1.3 |
        | :---: | :----: | :----: | :----: | :----: | :----: |
        | $\epsilon = 1$ | 184.7 | 78.24 | 66.1 | 55.3 | 46.3 |
        | $\epsilon = 0.1$ | 72.9 | 35.5 | 29.7 | 25.0 | 20.7 |
        | $\epsilon = 0.01$ | 6.77 | 3.98 | 3.35 | 2.84 | 2.34 |
        | $\epsilon = 0.0001$ | 1.40 | 1.30 | 1.69 | 不收敛 | 不收敛 |

        | 偏差最大值 | Jacobi | G-S | SOR w=1.1 | w=1.2 | w=1.3 |
        | :---: | :----: | :----: | :----: | :----: | :----: |
        | $\epsilon = 1$ | $0.61 \times 10^{-3}$ | $1.28 \times 10^{-3}$ | $1.1 \times 10^{-3}$ | $0.96 \times 10^{-3}$ | $0.83 \times 10^{-3}$ |
        | $\epsilon = 0.1$ | $0.009 | $0.009$ | $0.009$ | $0.009$ | $0.009$ |
        | $\epsilon = 0.01$ | $0.066$ | $0.066$ | $0.066$ | $0.066$ | $0.066$ |
        | $\epsilon = 0.0001$ | $0.005$ | $0.005$ | $0.005$ | 不收敛 | 不收敛 |
        简单分析一下，当 $\epsilon$ 接近 $1$ 时，收敛速度非常慢；
        当 $\epsilon$ 非常小时，对松弛因子的要求更高（因为对角元几乎不能占优）；
        当 $\epsilon$ 较小时，解析解会在 $x=0$ 附近发生“跳变”，使得离散化方程组的解近似为
        $$ y = (1-a) + ax $$
        表中的数据主要由跳变处（$i=1$）的偏差构成。提高跳变处附近的准确性，会成为优化的关键方向，不过这不属于解线性方程组的范畴
### 二、二阶偏微分方程的差分（拉普拉斯方程）
* 问题描述：
偏微分方程
    $$ 
    \begin{cases}
    & -\Delta u + g(x,y) u = f(x,y), & (x,y) \in D \\
    & u(x,y) = 1, & (x,y) \in \partial D\\
    & D = [0,1] \times [0,1]
    \end{cases}
    $$
同上题进行离散化后，同样得到大型稀疏方程组
    $$
    -u_{i-1,j} -u_{i,j-1} +(4+h^2 g(ih,jh))u_{i,j} -u_{i+1,j} -u_{i,j+1} = h^2 f(ih,jh)
    $$ $$1 \leq i,j \leq N-1$$
取$g(x,y) = e^{xy}, f(x,y) = x+y$进行计算
* 算法原理：
考虑相应的迭代格式，将其还原成系数矩阵，不难得到适应的迭代方程（源代码详见```Function.cpp```后半部分）：
    1. Jacobi
    $$
    -u_{i-1,j}^{(k)} -u_{i,j-1}^{(k)} +(4+h^2 g(ih,jh))u_{i,j}^{(k+1)} -u_{i+1,j}^{(k)} -u_{i,j+1}^{(k)} = h^2 f(ih,jh)
    $$
    2. G-S
    $$
    -u_{i-1,j}^{(k+1)} -u_{i,j-1}^{(k+1)} +(4+h^2 g(ih,jh))u_{i,j}^{(k+1)} -u_{i+1,j}^{(k)} -u_{i,j+1}^{(k)} = h^2 f(ih,jh)
    $$
    3. SOR
    $$
    -\omega (u_{i-1,j}^{(k+1)} +u_{i,j-1}^{(k+1)}) +(4+h^2 g(ih,jh))(u_{i,j}^{(k+1)}- (1-\omega)u_{i,j}^{(k)}) -\omega (u_{i+1,j}^{(k)} +u_{i,j+1}^{(k)}) = \omega h^2 f(ih,jh)
    $$
* 运行结果及分析：
    1. 对于不同的 $N$ ，解矩阵收敛的情况如下（最优松弛因子已验证，斜杠为未测结果，不代表算法不收敛）：

        | 迭代次数 | Jacobi | G-S | SOR w=1.25 | w=1.1 |
        | :---: | :----: | :----: | :----: | :----: |
        | $N=20$ | 54 | 31 | 18 | / | 
        | $N=40$ | 59 | 34 | 20 | 18 |
        | $N=60$ | 62 | 35 | 21 | 19 |

        | 耗时（单位：s） | Jacobi | G-S | SOR w=1.25 | w=1.3 |
        | :---: | :----: | :----: | :----: | :----: |
        | $N=20$ | 0.15 | 0.087 | 0.058 | / | 
        | $N=40$ | 0.57 | 0.34 | 0.21 | 0.19 |
        | $N=60$ | 1.29 | 0.72 | 0.49 | 0.43 |
    
    2. 由于解的规模比较庞大，不在报告中展示。求解该方程组的算法很快，如需查看可以直接运行程序

    相比前一题，理论上问题的规模提高了至少一个数量级，但是迭代法的求解反而更快。这主要得益于两方面的优势：针对系数矩阵的结构设计迭代格式，使得平均每次迭代所需时间大幅下降；问题的特殊性导致迭代次数很小（这里还没想清楚是为什么）。

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({ tex2jax: {inlineMath: [['$', '$']]}, messageStyle: "none" });
</script>