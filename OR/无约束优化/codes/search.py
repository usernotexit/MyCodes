# 给定函数f，出发点x和搜索方向d，进行非精确一维搜索，返回搜索步长
# input: f(x), x_0, d
# output: alpha
# parameters: rho in (0,1/2), sigma in (rho, 1)

# 参考文档：https://blog.csdn.net/weixin_43761124/article/details/107436454
# 添加了条件：若一维搜索 phi0 与 phi1 很接近而未完成搜索，则发出警告
# 此外还有一些优化空间：计算目标函数维度、梯度的解析式和Hesse阵的解析式重复次数很多，
#   这些重复计算可以省去，newton.py 同理

import numpy as np
from sympy import *
import re

dim_max = 10 # 目标函数最大维度
'''一些全局参数'''

class Line_Searcher:
    def __init__(self, rho=0.3, sigma=0.6):
        err_msg = 'please make sure rho in (0,1/2) and sigma in (rho, 1)'
        assert (rho>0 and rho<0.5 and sigma>rho and sigma<1) ,err_msg

        self.rho = rho
        self.sigma =sigma
        def fun(x):         # 为了避免忘记加函数出错，设置初始函数
            return 0*x[0]+0*x[1]+0*x[2]
        self.f = fun
        self.dim = 1

        self.failed = False # 为了避免反复输出 搜索失败 的信息，添加辅助元素

    def add_fun(self, fun):
        '''添加/更改目标函数，同时更新维度数据'''
        self.f = fun
        
        x_ = []
        for a in range(dim_max):
            x_.append(symbols('x'+str(a),real=True))  #设置符号变量
        f_ = self.f(x_)   # f的正则表达式f_
        dimensions = _dim(f_)
        self.dim = dimensions

    def WolfePowell(self, x, d, alpha_=1.):
        '''Wolfe-Powell非精确线性搜索，返回函数f在x处，方向d时的步长alpha'''
        maxk = 1000  #迭代上限
        k = 0
        phi0 = float(self.f(x))

        grad = _grad(self.f, x)
        phi0_ = float(np.dot(grad, d))
        #print(type(phi0_))
        
        a1 = 0
        a2 = alpha_
        alpha = (a1+a2)/2
        phi1 = phi0
        phi1_ = phi0_

        k = 0
        for k in range(maxk):    #限制迭代上限,避免时间太长
            phi = float(self.f(x + alpha * d))
            if phi <= phi1 + self.rho * alpha * phi1_:
                newx = x + alpha * d
                newgrad = _grad(self.f, newx)

                phi_ = np.dot(newgrad, d)
                
                if phi_ >= self.sigma*phi0_:
                    break
                else:
                    if abs(phi1_-phi_)<1e-10 and (not self.failed):
                        self.failed = True
                        print("Warning! One line search might have failed.",end='')
                        print("You can loop rho and sigma and retry.")
                        break                    
                    alpha_new = alpha + (alpha-a1) *phi_ / (phi1_-phi_)
                    a1 = alpha
                    alpha = alpha_new
                    phi1 = phi
                    phi1_ = phi_
            else:
                if abs((phi1-phi)-(a1-alpha)*phi1_)<1e-10 and (not self.failed):
                    self.failed = True
                    print("Warning! One line search might have failed.",end='')
                    print("You can loop rho and sigma and retry.")                    
                    break
                alpha_new = a1 + 0.5*(a1-alpha)**2*phi1_/((phi1-phi)-(a1-alpha)*phi1_)
                a2 = alpha
                alpha = alpha_new
        k = k + 1
        return alpha

def _dim(f_):
    '''利用正则表达式统计目标函数维度，输入为函数f的正则表达式f_'''
    dimension_set = []
    dimension_set = re.findall(r'x[0-9]\d*',str(f_))
    dimensions = len(set(dimension_set))
    return dimensions


def _grad(f, x):
    '''
    求f在x处的梯度（解析法），其步骤可以改进为 求函数维数+求解梯度解析式+解析式求值 ，
    减少求梯度解析式的次数
    '''
    x_ = []
    for a in range(dim_max):
        x_.append(symbols('x'+str(a),real=True))  #设置符号变量
    f_ = f(x_)   # f的正则表达式f_
    dimensions = _dim(f_)

    dx = []
    grad = []
    item={}
    for i in range(dimensions):
        item[x_[i]] = x[i]

    for a in range(dimensions):
        dx.append(diff(f_, symbols('x'+str(a),real=True))) #求偏导

        grad.append(dx[a].subs(item)) #求梯度 
    return grad