# 作业报告八
 - 陈泽高 PB20000302
## 代码结构：
* 本次作业含有三组代码，名称和功能如下
    * Homework8.cpp：运行，直接展示作业结果。
    * SVD.hpp：编写了 ```SVD``` 类用于实现奇异值分解。
    * Function.cpp/hpp： 包含 givens变换、矩阵减法、乘法等功能函数，以及通用的 ```output()``` 输出函数，可以指定输出的数据类型（数组、二维数组）以及精度。

    此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码；在 ```Homework8.cpp``` 中可以调整部分参数，查看需要检查的解

## 奇异值分解（SVD）
 * 问题描述：略
 * 算法原理：
    1. 对于矩阵 $A \in \R^{m\times n}$ ，存在正交阵 $P \in \R^{m\times m}$, $Q \in \R^{n\times n}$ 以及对角阵 $\Sigma \in \R^{m\times n} $ 满足 $$ A = P \Sigma Q$$ 其中 $\Sigma$ 的对角元称为 $A$ 的奇异值。在不考虑奇异值顺序的情况下，奇异值分解具有唯一性

    2. 对 $A^{T} A$ 隐式地应用对称QR方法：对 $A$ 二对角化，即计算两个正交矩阵 $U$ 和 $V$ 使得 $$ UAV= \begin{bmatrix}
    B \\
    0 \\
    \end{bmatrix} ,
    B = \begin{bmatrix}
    \alpha_{1} & \beta_{1} & \\
     & \ddots & \ddots & \\
     & & \ddots & \beta_{n-1} \\
     & & & \alpha_{n}
    \end{bmatrix} $$ 这个过程相当于将 $A^{T} A$ 三对角化，具体算法见教材 ```算法7.6.1``` ，代码见SVD.hpp中的函数 ``` _hessenberg(&P, &S, &Q)``` 。再对 $B^{T}B$ 应用隐式QR迭代，得到$$ B^{T}B = (Q^{T} \Sigma^{T} P^{T})(P \Sigma Q) $$ 其中 $ \Sigma^{T}\Sigma $ 是对角阵，即 $\Sigma$ 是对角阵，具体算法见教材 ```算法7.6.2``` ，代码见SVD.hpp中的函数 ```_qr(&P,&S,&Q,r,l)``` 以及 ```_adapt(&S, &P, r, l, pos)```， 至此完成对A的奇异值分解。
* 实验结果及分析：
    待分解矩阵 $A$ 在文件 ```svddata.txt``` 中，奇异值如下：
    <center>
    <img src="screenshots/result.jpg" style="zoom: 80%;" />
    </center>
    
    迭代次数为 $17$ 在合理范围内，奇异值序列的精度经过验证，$P$,$Q$ 作为正交阵的误差以及 $P\Sigma Q-A$ 的最大元素如图。总之，SVD分解精度符合预期。更具体的输出结果详见附件程序。