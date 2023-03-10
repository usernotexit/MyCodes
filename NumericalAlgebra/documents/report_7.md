# 作业报告七
 - 陈泽高 PB20000302
## 代码结构：
* 本次作业含有四组代码，名称和功能如下
    * Homework7.cpp：运行，直接展示作业结果
    * Symmetric.hpp：编写了 ```Threshold_Jacobi``` 类用于实现过关Jacobi方法求特征值与特征向量，以及 ```divided``` 类用于实现二分法求实对称三对角阵的指定特征值，为了方便调用，笔者暂时将**反幂法**内嵌到 ```divided``` 类中，可以根据未来的需求将其单独打包。
    * Solve.cpp/hpp： 实现了 ```Solver``` 类，这是以前作业使用的类，本次作业在实现反幂法时利用它来对系数矩阵做预分解，提高运算速度。本次使用列主元高斯分解，如有要求仅查看相应部分即可。
    * Function.cpp/hpp： 包含 ```Solver``` 类的功能函数，以及通用的 ```output()``` 输出函数，可以指定输出的数据类型（数组、二维数组）以及精度。

    此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码；在 ```Homework7.cpp``` 中可以调整部分参数，查看需要检查的解

## 一、过关Jacobi法求解实对称阵的正交对角化
 * 问题描述：略
 * 算法原理：
    1. 用非对角“范数”衡量矩阵A到对角阵的接近程度：
    $$
        E(A)=(\sum_{i=1}^{n} \sum_{j=1,j\ne i} \alpha_{i j}^{2})
    $$
    Jacobi方法通过Givens变换将$E(A)$逐渐降低为$0$（详见教材及课件），这样矩阵$A$将会收敛于对角阵$D$，而由于Givens变换的正交性，$D$的对角元就是$A$的特征值；相应的，将Givens变换逐次累乘，得到的矩阵$Q$的每一列都是相应特征值的特征向量
    
    2. 我们希望每次都找到绝对值最大的非对角元$\alpha_{i j}$，这样能使得收敛速度最大化，但是每次搜索的复杂度是平方级（每次Givens变换的复杂度仅为线性），时间开销过大。为此，我们采用过关Jacobi法：先确定一个阈值$\delta$，每次搜索对大于该阈值的非对角元$\alpha_{i j}$所在行、列进行Givens变换，当所有非对角元的绝对值都小于该阈值后，降低阈值并重复上述过程直至阈值足够小。本次作业阈值的选取如下：
    $$
        \delta_{0}=E(A) \\
        \delta_{n}=\delta_{n-1}/n\\
    $$
    
    3. Jacobi方法的显著优点是并行效率高，因此开启VS的```release```模式可以明显提高运算速度

 * 运行结果及分析：
    1. $n=50$
    由于保留了非零的非对角元，因此截图可能显示不全，详见附件程序的运行结果
    $A_k$的一部分：
    <center>
	<img src="screenshots/exe1_A_50.jpg" style="zoom: 40%;" />
    </center>
    
    $Q$ 的截图
    <center>
    <img src="screenshots/exe1_Q_50.jpg" style="zoom: 40%;" />
    <img src="screenshots/exe1_Q2_50.jpg" style="zoom: 40%;" />
    <img src="screenshots/exe1_Q3_50.jpg" style="zoom: 40%;" />
    <img src="screenshots/exe1_Q4_50.jpg" style="zoom: 40%;" />
    <img src="screenshots/exe1_Q5_50.jpg" style="zoom: 40%;" />
    </center>

    2. $n=50,60,70,80,90,100$ 的特征值全体（按顺序展示）
    <center>
    <img src="screenshots/exe1_50.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe1_60.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe1_70.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe1_80.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe1_90.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe1_100.jpg" style="zoom: 80%;" />
    </center>

    笔者用Mathematica软件对 $n=50,60$ 的情况进行了验证，结果的准确度较高
    <center>
    <img src="screenshots/exe1_50_refer.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe1_60_refer.jpg" style="zoom: 80%;" />
    </center>

## 二、二分法和反幂法
 * 问题描述：用二分法求实对称三对角阵指定特征值，再用反幂法求解该特征值对应的特征向量
 * 算法原理：
    1. 对于矩阵$T$，记 $p_{i}(\lambda)$ 为 $T-\lambda I$ 的$i$阶顺序主子式，满足如下递推公式：
    $$
        p_{0}(\lambda)=1,  p_{1}(\lambda)=\alpha_1 - \lambda \\
        p_{i}(\lambda)=(\alpha_i - \lambda)p_{i-1}(\lambda) - \beta_{i}^{2} p_{i-2}(\lambda) \\
        i=2,3,\cdots, n \\
    $$
    由于 $T$ 的特征值均为实数，$p_{i}(\lambda)$ 的根都是实的。对于不可约对称三对角阵，通过计算 $p_{i}(\lambda)$ 的变号数（$i=0,1,\cdots,n$）可以得到 $T$ 在区间 $(-\infty, \lambda]$ 中的特征值个数，结合二分法可以快速将目标特征值确定在一个极小的范围内，从而获得近似解。二分法的初始端点选为 $\parallel T\parallel_{\infty}$ 和 $ -\parallel T\parallel_{\infty}$

    2. 对于矩阵$T$及已知的特征值$\lambda$，可以对 $(T-\lambda I)^{-1}$ 应用幂法，一般来说只需迭代一步，得到的向量就与 $\lambda$ 对应的特征向量非常接近。本次作业笔者采用列主元Gauss分解的方法求解线性方程组 $(T-\lambda I)u=v$ ，其中$v$为随机的初始向量，得到$u$之后将其归一化，再替代$v$，重复前述过程若干次，最终得到归一化的特征向量$u$

 * 结果展示与分析：
    <center>
    <img src="screenshots/exe2_min.jpg" style="zoom: 80%;" />
    <img src="screenshots/exe2_max.jpg" style="zoom: 80%;" />
    </center>