import numpy as np
class LogisticRegression:
    
    def __init__(self, penalty="l2", fit_intercept=True, la=0.0):
        # gamma（好像是代码中self.la的作用）和penalty（只写了l2正则项）
        err_msg = "penalty must be 'l1' or 'l2', but got: {}".format(penalty)
        assert penalty in ["l2", "l1"], err_msg

        self.w = np.random.random((11,1)) + 4   # 初始值完全中立，则初始loss就很小，看不出明显提升
        self.w_hat = 10 * self.w
        self.fit_intercept = fit_intercept
        self.la = la
        self.penalty = penalty

        if not fit_intercept:
            self.b = 0.0
            self.w_hat = np.row_stack((self.w, self.b))
        
        self.losscurve = np.array([0,0])
        self.loss = np.array(0)

    def sigmoid(self, x):
        """The logistic sigmoid function"""
        return 1.0/(1 + np.exp(-x.transpose() * self.w_hat))# self.w_hat

    def draw_loss(self, iter):
        #epoch = (iter-1)/100        
        loss_add = np.array([iter, self.loss])
        self.losscurve = np.column_stack((self.losscurve, loss_add))

    def Reshape(self, X, add=True):
        '''重构数据集X,y等以便进行矩阵操作'''
        X = np.mat(X).transpose()
        if add & (not self.fit_intercept):
            Xadd = np.ones(X.shape[1])      # bias
            X = np.row_stack((X, Xadd))     # 得到x_hat
        return X

    def fit(self, X, y, lr=0.001, max_iter=100000, tol=1e-5, draw_losscurve=False):
        """
        Fit the regression coefficients via gradient descent or other methods 
        lr: 设置移动步长(学习率)
        tol: tolerance 停止标准，误差不超过tol时，停止进一步的计算
        max_iter: 最大循环次数，当达到该次数时，即使偏差没有小于tol，仍然停止循环
        """        
        X = self.Reshape(X)
        y = self.Reshape(y, add=False)
        num = y.shape[0]
        self.max_loss = 0.0 # 局部最大loss

        for iter in range(1,max_iter):
            '''开始循环，梯度下降法'''
            h = self.sigmoid(X)         # X即为数据集x_hat
            grad = y - h

            if self.penalty=='l2':
                self.w_hat += lr * X * grad - lr * self.la * self.w_hat
            else:
                self.w_hat += lr * X * grad

            self.loss = (- y.transpose() * np.log(h + 1e-5) - (1-y.transpose()) * np.log(1-h + 1e-5))/num
            # print(self.loss) #测试用

            if draw_losscurve :#& (iter%10 == 1):
                self.draw_loss(iter)
            
            '''停止条件'''
            if self.loss < tol:
                print('break1')
                break
            if self.loss > self.max_loss:
                self.max_loss = self.loss
            if iter%1000 == 0:
                '''一轮迭代100次'''
                if (self.max_loss - self.loss) < tol/10:
                    '''缩短步长以获得更准确结果 + early stop'''#偷个懒用tol做条件
                    print('break2')
                    break
                    #lr = lr/2
                    # print(lr)
                    #if lr < 1e-8:   #最小步长，这个也可以写成lr_max/100，或者作为参数传入
                    #    print('break2')
                    #    break
                self.max_loss = 0
        #print("Training over")
        
    def predict(self, X, y):
        """
        Use the trained model to generate prediction probabilities on a new
        collection of data points.
        """
        X = self.Reshape(X)
        y = self.Reshape(y, add=False)
        
        h = np.array(self.sigmoid(X))
        h[h>0.5] = 1.0
        h[h<=0.5] = 0.0
        accu = 1.0 - ((h-y).transpose() * (h-y))/(y.shape[0])
        return accu
    
    def predict_1(self, X):
        X = self.Reshape(X)
        
        h = np.array(self.sigmoid(X))
        h[h>0.5] = 1.0
        h[h<=0.5] = 0.0
        return h