# 作业报告一
## 陈泽高 PB20000302
### 代码结构：
* 本次作业含有四组代码，名称和功能如下
    * Homework1.cpp：运行终端、用户界面，展示作业结果
    * Exercise.cpp/hpp：其中每个函数各自完成一个题目，构造方程并调用Solve.h（库）的函数求解方程
    * Solve.cpp/hpp：调用Function.h的线性方程组求解器，给出求解时间和误差
    * Function.cpp/hpp：Gauss消元法、Cholesky分解法等涉及矩阵运算的函数、还有检查函数check()
	
	此外，在注释中还保留了若干测试痕迹，审阅时可以利用它们来测试代码
### 一、Gauss消去法
* 问题描述：将不选主元、全主元、列主元的Gauss消去法编写成通用子程序，求解84阶方程组
并将结果与精确解进行比较
* 算法原理：三种Gauss消去法，前代法和回代法解三角形方程组，详见教材描述。
Gauss消去法核心代码如下：（着重突出算法，假设）
    1. 不选主元
    ```
    // （不选主元）高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
    // 参考教材上的算法
	for(k=0; k<N; k++){
	    	// 计算L的第(k+1)列 
		for(j=k+1; j<N; j++){
    		A[j][k] = A[j][k] / A[k][k];
	   	}

		// 让剩下部分适应前面的变换，自动计算U的第(k+1)行 
		for(i=k+1; i<N; i++){
		    for(j=k+1; j<N; j++){
			    A[i][j] = A[i][j] - A[i][k]*A[k][j];
		    }
		}
	}
    ```
    2. 全主元
    ```
    // 全主元高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
	// 参考教材上的算法，PAQ = LU  
    for(k=0; k<N; k++){
		// 选取主元 
		p=k, q=k;		//记录主元位置 
		for(i=k; i<N; i++){
			for(j=k; j<N; j++){
				if(|A[i][j]| > |A[p][q]|) {
					p = i;
					q = j;
				}
			}
		}
		for(i=0; i<N; i++){	// exchange row k with row p 
			swamp = A[k][i]; A[k][i] = A[p][i]; A[p][i] = swamp;
		}
		for(j=0; j<N; j++){	//exchange column k with column q
			swamp = A[k][i]; A[k][i] = A[p][i]; A[p][i] = swamp;
		}
		u[k] = p; v[k] = q;//记录矩阵P,Q
		
		// 接下来和普通的高斯消元法一样 
		for(j=k+1; j<N; j++){
			A[j][k] = A[j][k] / A[k][k];
		}
		for(i=k+1; i<N; i++){
			for(j=k+1; j<N; j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
			}
		}
	}
    ```
    3. 列主元
    ```
    // 列主元高斯消去法：对A进行LU分解，并分别存储在A的上半部分（含对角线）和下半部分（不含对角线） 
	// 参考教材上的算法，PA = LU  
    for(k=0; k<N; k++){
		// 选取主元 
		p=k, q=k;		//记录主元位置 
		for(i=k; i<N; i++){
			if(|A[i][j]| > |A[p][q]|) {
				p = i;
			}
		}
		for(i=0; i<N; i++){	// exchange row k with row p 
			swamp = A[k][i]; A[k][i] = A[p][i]; A[p][i] = swamp;
		}
		u[k] = p;	记录矩阵P
		
		// 接下来和普通的高斯消元法一样 
		for(j=k+1; j<N; j++){
			A[j][k] = A[j][k] / A[k][k];
		}
		for(i=k+1; i<N; i++){
			for(j=k+1; j<N; j++){
				A[i][j] = A[i][j] - A[i][k]*A[k][j];
			}
		}
	}
    ```

* 运行结果与对比：解84阶方程组
    <center>

	***equations***

	<img src="images_and_screenshots/Equations1.jpg" style="zoom: 70%;" />
	
	***results***
    
	<img src="images_and_screenshots/1_gauss_1.jpg"  style="zoom: 70%;" />
    <img src="images_and_screenshots/1_gauss_2.jpg" style="zoom: 70%;" />
    <img src="images_and_screenshots/1_gauss_3.jpg" style="zoom: 70%;" />	

	</center>
	
	可见不选主元运算时间最短（显然），但准确度在第40个分量后便大幅下降，且偏差$err$较大；全主元和列主元方法解的质量相近，但列主元方法省下了许多时间


### 二、Cholesky分解法
* 问题描述：将平方根法和改进的平方根法编写成通用的子程序，求解正定对称阵的线性方程组
* 算法描述：正定对称阵的平方根分解法和改进的平方根分解法，详见教材描述。
为了便于调用前代法和回代法，还将只下三角的部分拷贝到上三角，详见原代码Function.cpp
两种方法的核心代码如下：
    1. 平方根法
    ```
	// 对正定对称阵用cholesky算法分解，把结果储存在矩阵的下三角部分 
	for(k=0; k<N; k++){
		A[k][k] = sqrt(A[k][k]);
		for(j=k+1; j<N; j++){
			A[j][k] = A[j][k]/A[k][k];
		}
		for(j=k+1; j<N; j++){
			for(i=j; i<N; i++){
				A[i][j] = A[i][j] - A[i][k]*A[j][k];
			}
		}
	}
    ```

    2. 改进的平方根法
    ```
	// 对正定对称阵用改进的cholesky算法分解，把结果储存在矩阵的下三角部分 
	for(j=0; j<N; j++){
		for(i=0; i<j; i++){
			v[i] = A[j][i]*A[i][i];
			A[j][j] = A[j][j] - A[j][i]*v[i];
		}
		for(i=j+1; i<N; i++){
			for(k=0; k<j; k++){
				A[i][j] = A[i][j] - A[i][k]*v[k];
			}
			A[i][j] = A[i][j]/A[j][j];
		}
	}
    ```
* 运行结果与对比：
	* 解100阶方程组，向量$b$随机选取
	本次实验中，选出的$b_i$是$[0,100]$内的随机双精度浮点数
	<center>

	***Matrix***
	
	<img src="images_and_screenshots/Equations2.jpg" style="zoom: 70%;" />
	
	***results***
    
	<img src="images_and_screenshots/2_cholesky_1.jpg"  style="zoom: 70%;" />
    <img src="images_and_screenshots/2_cholesky_2.jpg" style="zoom: 70%;" />

	</center>

	可见修改后的平方根法无论在精度还是速度方面都胜过原平方根法————开方的运算慢，且造成的偏差较大

	* 解40阶Hilbert矩阵对应的方程组，将精确解设定为1（所有分量都是1）
	<center>

	***results***
    
	<img src="images_and_screenshots/2_cholesky_3.jpg"  style="zoom: 70%;" />
	</center>

	在解Hilbert矩阵的方程时，当阶数较大时，两种方法都不可避免地出现了偏差（尽管$err$还很小）；当阶数增长到一定程度（本机是13），平方根法便会超出double的取值范围。
### 三、Gauss消元法与Cholesky分解法对比
* 问题描述：用Gauss消元法和平方根法分别求解对称正定方程组，并对比求解时间和精度等
* 算法描述：
    1. 两种方法的代码见前两个板块；
    2. 精度估计：对于方程$Ax=b$的机器解$x$和精确解$x_0$，一种可以参考的指标是：<center> $err=|Ax-b|_2$ 
        
    $err$越小，$|x-x_0|_2$越小，即得到的解与精确解更接近。

* 运行结果与分析：
	* 解100阶方程组，向量$b$随机选取
	本次实验中，选出的$b_i$是$[0,100]$内的随机双精度浮点数
	<center>

	***Matrix***
	
	<img src="images_and_screenshots/Equations2.jpg" style="zoom: 70%;" />
	
	***results***
    
	<img src="images_and_screenshots/3_cholesky_rand_1.jpg"  style="zoom: 70%;" />
    <img src="images_and_screenshots/3_cholesky_rand_2.jpg"  style="zoom: 70%;" />
    <img src="images_and_screenshots/3_gauss_rand_1.jpg"  style="zoom: 70%;" />
	<img src="images_and_screenshots/3_gauss_rand_2.jpg"  style="zoom: 70%;" />
    <img src="images_and_screenshots/3_gauss_rand_3.jpg"  style="zoom: 70%;" />
	</center>

	1. 求解对角占优阵的方程时，五种方法都有良好的表现；
	2. 平方根法以及其改进方法都远快于高斯消元法，精度方面则相差不大。
	3. 随着矩阵阶数的提升，平方根法将会有更好的表现


	* 解40阶Hilbert矩阵对应的方程组，将精确解设定为1（所有分量都是1）
	<center>

	***results***

    </center>
	
	<center>
		<img src="images_and_screenshots/3_cholesky_Hilbert.jpg"  style="zoom: 70%;" />
		<img src="images_and_screenshots/3_gauss_Hilbert.jpg"  style="zoom: 70%;" />
	</center>

	不难看出，现有的五种方法在求解Hilbert矩阵的方程组时，都出现了很大的偏离（尽管$err$不大）。这也说明该矩阵在向量的2范数下具有某种“不稳定性”
	