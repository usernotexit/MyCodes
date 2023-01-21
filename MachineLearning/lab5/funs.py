import time as t
import numpy as np
import pandas as pd

class SVM:
    '''
    支持向量机模型（二分类）,采用SGD或者SMO算法求解\n
    建议调用SGD，因为从之前的测试来看该算法更优
    '''
    def __init__(self, lr=1e-4, tol=1e-4, C=10, max_iter=10):
        """
        lr:学习率
        tol:当步长小于该值，训练停止
        C:正则项的倒数，C越小，KKT条件判断越严格
        max_iter:
        """
        #self.dim = dim
        self.lr = lr
        self.tol = tol
        self.C = C
        self.max_iter = max_iter

    def fit(self, X, y):
        """
        Fit the coefficients via SGD (mini-batch) method
        """
        self.dim = X.shape[1]
        
        w = np.zeros((1, self.dim))
        b = 0
        #cost = []

        start = t.time()
        for j in range(self.max_iter):
            arr=np.array([i for i in range(X.shape[0])])
            np.random.shuffle(arr)
            for i in arr:
                if y[i] * (X[i].dot(w.T)) >= 1:
                    gradient = 0
                    grad_b = 0
                else :
                    gradient = -y[i] * X[i]
                    grad_b = -y[i]
                gradient = w + (self.C * gradient)
                grad_b = self.C * grad_b
                w = w - self.lr * gradient
                b = b - self.lr * grad_b
                #cost.append(self.cost(X, y, w, b))

            if np.linalg.norm(self.lr * gradient) < self.tol:
                break
     
        self.w = w
        self.b = b

        end = t.time()
        self.train_time = end - start
        return #cost

    def predict(self, X):
        """
        Use the trained model to generate prediction on a
        collection of data points.\n
        X is n*dim 's array
        """
        y_hat = np.dot(X, self.w.T) + self.b
        y_hat[y_hat>=0] = 1.0#
        y_hat[y_hat<-0] = 0.0#
        return y_hat

    def predict_proba(self, X):
        """
        Use the trained model to generate prediction probabilities on a
        collection of data points.\n
        X is n*dim 's array
        """
        y_hat = np.dot(X, self.w.T) + self.b
        y_hat = np.array(self.sigmoid(y_hat))
        return y_hat
    #def decision_function(self,X):

    def score(self, X, y):
        """
        accuracy score
        """
        y_hat = self.predict(X)
        miss = y-y_hat
        return 1 - miss.T@miss / len(y)

    def sigmoid(self, x):
        """The logistic sigmoid function"""
        return 1.0/(1 + np.exp(-x))# self.w_hat

    def cost_func(self, X, y, w, b):
        """The cost function"""
        result = 0
        for i in range(X.shape[0]):
             y_pre = X[i].dot(w.T) + b
             if y[i] * y_pre >= 1 :
                 continue
             else :
                 result = result + (1 - y[i] * y_pre)
        result = result / X.shape[0]
        result = self.C * result + (np.linalg.norm(w) ** 2) / 2

        return result
    
    def get_params(self, deep=False):
        '''获取模型参数'''
        dic = dict()
        dic['lr'] = self.lr
        dic['tol'] = self.tol
        dic['max_iter'] = self.max_iter
        dic['C'] = self.C
        if deep:
            pass
        return dic

class LogisticRegression:
    '''逻辑回归模型'''
    def __init__(self, penalty="l2", fit_intercept=True, la=0.0, lr=0.001, max_iter=100000, tol=1e-5):
        '''
        fit_intercept：是否固定截距
        lr: 设置移动步长(学习率)
        tol: tolerance 停止标准，误差不超过tol时，停止进一步的计算
        max_iter: 最大循环次数，当达到该次数时，即使偏差没有小于tol，仍然停止循环
        la：正则项参数
        '''
        # gamma（好像是代码中self.la的作用）和penalty（只写了l2正则项）
        err_msg = "penalty must be 'l1' or 'l2', but got: {}".format(penalty)
        assert penalty in ["l2", "l1"], err_msg

        self.fit_intercept = fit_intercept
        self.la = la
        self.penalty = penalty
        self.lr = lr
        self.max_iter = max_iter
        self.tol = tol
        
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

    def fit(self, X, y, draw_losscurve=False):
        """
        Fit the regression coefficients via gradient descent or other methods 
        """
        n, m = X.shape
        self.w = np.random.random((m,1)) + 4   # 初始值完全中立，则初始loss就很小，看不出明显提升
        self.w_hat = 10 * self.w
        if not self.fit_intercept:
            self.b = 0.0
            self.w_hat = np.row_stack((self.w, self.b))
        
        X = self.Reshape(X)
        y = self.Reshape(y, add=False)
        num = y.shape[0]
        self.max_loss = 0.0 # 局部最大loss

        for iter in range(1,self.max_iter):
            '''开始循环，梯度下降法'''
            h = self.sigmoid(X)         # X即为数据集x_hat
            grad = y - h

            if self.penalty=='l2':
                self.w_hat += self.lr * X * grad - self.lr * self.la * self.w_hat
            else:
                self.w_hat += self.lr * X * grad

            self.loss = (- y.transpose() * np.log(h + 1e-5) - (1-y.transpose()) * np.log(1-h + 1e-5))/num
            # print(self.loss) #测试用

            if draw_losscurve :#& (iter%10 == 1):
                self.draw_loss(iter)
            
            '''停止条件'''
            if self.loss < self.tol:
                print('break1')
                break
            if self.loss > self.max_loss:
                self.max_loss = self.loss
            if iter%1000 == 0:
                '''一轮迭代100次'''
                if (self.max_loss - self.loss) < self.tol/10:
                    '''缩短步长以获得更准确结果 + early stop'''#偷个懒用tol做条件
                    print('break2')
                    break
                    #lr = lr/2
                    # print(lr)
                    #if lr < 1e-8:   #最小步长，这个也可以写成lr_max/100，或者作为参数传入
                    #    print('break2')
                    #    break
                self.max_loss = 0

    def predict(self, X):
        X = self.Reshape(X)
        
        h = np.array(self.sigmoid(X))
        h[h>0.5] = 1.0
        h[h<=0.5] = 0.0
        return h
        
    def predict_proba(self, X):
        X = self.Reshape(X)
        h = np.array(self.sigmoid(X))
        return h
        
    def score(self, X, y):
        """
        准确率得分，返回判定准确率
        """
        #X = self.Reshape(X)
        y = self.Reshape(y, add=False)
        
        h = self.predict(X)
        accu = 1.0 - ((h-y).transpose() * (h-y))/(y.shape[0])
        return accu
    
    def get_params(self, deep=False):
        '''获取模型参数'''
        dic = dict()
        # dic['w'] = self.w_hat
        dic['fit_intercept'] = self.fit_intercept
        dic['la'] = self.la
        dic['lr'] = self.lr
        dic['tol'] = self.tol
        dic['max_iter'] = self.max_iter
        dic['penalty'] = self.penalty
        if deep:
            pass
        return dic


def distanceNorm(Norm,D_value):
    '''计算向量范数'''    
    # Norm for distance
    if Norm == '1':
        counter = np.absolute(D_value)
        counter = np.sum(counter)
    elif Norm == '2':
        counter = np.power(D_value,2)
        counter = np.sum(counter)
        counter = np.sqrt(counter)
    elif Norm == 'Infinity':
        counter = np.absolute(D_value)
        counter = np.max(counter)
    else:
        raise Exception('We will program this later......')
    return counter

#原文链接：https://blog.csdn.net/kryolith/article/details/40203817

class decision_tree:
    def __init__(self, depth=10):
        #self.dim = dim #只存在根节点即可
        self.maxdepth = depth
        self.is_leaf = True
    
    def fit(self, X, g, lamb=0., gamma=0., delta=0.):#生成决策树
        [m, dim] = X.shape
        self.dim = dim  #只存在根节点即可

        Obj, w = self.Obj(g, lamb, gamma)

        self.fit_(X, g, w, Obj, lamb, gamma, delta)

    def fit_(self, X, g, w, Obj, lamb, gamma, delta):#
        '''生成决策树的后续流程 g:残差g=y w:预估权重 Obj:单叶得分'''
        if self.maxdepth<=0 or g.size<5:
            self.w = w
            self.penal = gamma + 1/2 *lamb* w*w
            #print(self.penal)
            return

        '''获取最佳划分及其得分'''
        a, Na, gain, wL, wR, ObjL, ObjR = self.get_attri(X, g, w, Obj, lamb, gamma, delta)

        if gain<=delta: # delta>0
            self.w = w
            self.penal = gamma + 1/2 *lamb* w*w
            #print(self.penal)
            return
        else:
            self.is_leaf = False
            self.attri = a
            self.Na = Na #判定属性及其阈值
            XL, XR, gL, gR= self.split(a, Na, X, g)
            self.Lchild = decision_tree(self.maxdepth-1)
            self.Rchild = decision_tree(self.maxdepth-1)
            self.Lchild.fit_(XL, gL, wL, ObjL, lamb, gamma, delta)#wL
            self.Rchild.fit_(XR, gR, wR, ObjR, lamb, gamma, delta)#wR
        
    def get_attri(self, X, g, w, Obj, lamb, gamma, delta):
        '''try to get the optimal attribute and point to split'''
        # gain = Obj(P) - Obj(L) - Obj(R) # Obj(P)=sum(g)^2*(...)+gamma
        [m, dim] = X.shape
        attri = 0
        [gain_max, Nattri, wL, wR, ObjL, ObjR] = [delta, 0., 0., 0., 0., 0.]

        for a in range(dim):
            g = g[X[:,a].argsort()]#按属性排序
            X = X[X[:,a].argsort()]
            '''s:等分参数，用于时间优化，内置为10~100'''
            s = 10
            for split in range(s-1):
                '''为了减少计算量，改为十等分'''
                if m>s:
                    Na = (X[int(m/s)*(split+1), a] + X[int(m/s)*(split+1)+1, a])/2
                elif split<m-1:
                    Na = (X[split, a] + X[split+1, a]) / 2            
                else:
                    continue        
                XL, XR, gL, gR = self.split(a, Na, X, g)
                if gL.size<=0 or gR.size<=0:
                    #print('skip')
                    continue
                ObjL_tmp, wL_tmp = self.Obj(gL, lamb, gamma)
                ObjR_tmp, wR_tmp = self.Obj(gR, lamb, gamma)
                
                #mL = gL.size
                #mR = gR.size
                #gain = - sum(g.reshape(m))*sum(g.reshape(m))/ (m) + sum(gL.reshape(mL))*sum(gL.reshape(mL))/ (mL)\
                #    + sum(gR.reshape(mR))*sum(gR.reshape(mR))/ (mR)
                
                gain = Obj - ObjL_tmp - ObjR_tmp
                #if gain<0:
                #    print(gain) 
                #else:
                #    print('yes!')
             
                if gain>gain_max:
                    #print(gain)
                    attri = a
                    gain_max = gain
                    Nattri = Na            
                    wL = wL_tmp
                    wR = wR_tmp
                    ObjL = ObjL_tmp
                    ObjR = ObjR_tmp

        return attri, Nattri, gain_max, wL, wR, ObjL, ObjR

    def Obj(self, g, lamb, gamma): #实际上与结点无关
        '''return the Obj and optimal weight of the node'''
        m = g.size
        '''公式检查点'''
        w = sum(g.reshape(m)) / (m + lamb/2)
        #print(w)
        Obj = - sum(g.reshape(m))*sum(g.reshape(m))/ (m + lamb/2) + gamma
        return Obj, w

    def Get_penal(self):
        '''获取树的惩罚项:叶子结点*gamma + |w|^2求和'''
        penal = 0.
        if self.is_leaf:
            penal = self.penal
        else:
            penal = self.Lchild.Get_penal() + self.Rchild.Get_penal()
        return penal
        
    def split(self, a, Na, X, g):
        XL = X[X[:,a]<Na,:]
        XR = X[X[:,a]>=Na,:]
        gL = g[X[:,a]<Na,:]
        gR = g[X[:,a]>=Na,:]
        return XL, XR, gL, gR
        
    def predict(self, x):
        if self.is_leaf:
            w_hat = self.w
        elif x[self.attri]<self.Na:
            w_hat = self.Lchild.predict(x)
        else:
            w_hat = self.Rchild.predict(x)
        return w_hat

    def predict_(self, X):
        m = X.shape[0]
        # y_hat VS w_hat
        w_hat = []
        for i in range(m):
            w_hat.append(self.predict(x=X[i,:]))
        return w_hat
