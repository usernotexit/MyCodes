# 牛顿法求解无约束最优化问题，可拓展到其他方法
# input: f(x), x_0
# output: x_min = argmin_{x} f(x)
# parameters (for line search): rho in (0,1/2), sigma in (rho, 1)

import numpy as np
from sympy import *
import re
#import torch

from search import Line_Searcher
from search import _grad, _dim, dim_max

class opt:
    def __init__(self, rho=0.3, sigma=0.6):
        err_msg = 'please make sure rho in (0,1/2) and sigma in (rho, 1)'
        assert (rho>0 and rho<0.5 and sigma>rho and sigma<1) ,err_msg

        self.rho = rho
        self.sigma =sigma
        def fun(x):         # 为了避免忘记加函数出错，设置初始函数
            return 0*x[0]+0*x[1]+0*x[2]
        self.f = fun
        self.dim = 1

    def add_fun(self, fun):
        '''添加/更改目标函数，同时更新维度数据'''
        self.f = fun
        
        x_ = []
        for a in range(dim_max):
            x_.append(symbols('x'+str(a),real=True))  #设置符号变量
        f_ = self.f(x_)   # f的正则表达式f_
        dimensions = _dim(f_)
        self.dim = dimensions

    def opt_newtown(self, x, max_iter=1000, grad_gap=1e-12):
        ''' 用牛顿法+非精确一维搜索完成无约束最优化 '''
        searcher = Line_Searcher(self.rho, self.sigma)
        searcher.add_fun(self.f)

        self.stop = 0  # 停机条件 0:达到最大迭代次数
        for iter in range(max_iter):    # 最大迭代次数，可自行设置
            '''搜索方向'''
            grad = _grad(self.f, x)
            grad = list(map(float, grad))
            grad = np.array(grad)
            if np.linalg.norm(grad) < grad_gap:    # 停机条件 1:梯度阻尼
                self.stop = 1
                break

            hesse = _hesse(self.f, x)
            hesse = [list(map(float, e)) for e in hesse]
            hesse = np.array(hesse)

            if abs(np.linalg.det(hesse)) < 1e-10*np.linalg.norm(hesse, np.inf):
                direc = -grad   # hesse可能较为病态，不采用牛顿步
                print("采用梯度法一次")
            else:
                direc = -np.linalg.inv(hesse) @ grad
            
            '''一维搜索'''
            print(type(direc[0]))####################################33
            alpha = searcher.WolfePowell(x, direc, 1.)
            x = x + alpha*direc
        
        if self.stop == 0:
            print("触发停机条件：达到最大迭代次数，注意此时可能未找到最优解，", end='')
            print("请增加迭代次数或增大梯度阈值重试。")
        elif self.stop == 1:
            print("触发停机条件：梯度小于阈值，运行正常")
        
        return x

    def opt_dfp(self, x0, max_iter=1000, delta=1e-10):
        ''' 用拟牛顿法+非精确一维搜索完成无约束最优化 '''
        searcher = Line_Searcher(self.rho, self.sigma)
        searcher.add_fun(self.f)

        x = x0
        d = np.array([-float(e) for e in _grad(self.f, x)])
        H = np.eye(len(x)) # H0 = I
        self.stop = 0  # 停机条件 0:达到最大迭代次数

        for iter in range(max_iter):
            alpha = searcher.WolfePowell(x, d, 1.)
            x_new = x + alpha * d
            
            if np.linalg.norm(_grad(self.f, x_new), np.inf) < delta: #np.linalg.norm(np.array([float(e) for e in x_new-x])) < delta:
                '''修改停机条件'''
                self.stop = 1 # 停机条件 1:梯度阻尼
                break
            s = (x_new - x).reshape(len(x),1)
            y = (_grad(self.f, x_new) - np.array([float(e) for e in _grad(self.f, x)])).reshape(len(x),1)
            #print(np.dot(y, y.T))
            H = H - np.dot(np.dot(H, np.dot(y, y.T)), H) / np.dot(np.dot(y.T, H), y) + np.dot(s, s.T) / np.dot(s.T, y)
            d = -np.dot(H, _grad(self.f, x_new)) # 搜索方向
            x = x_new
        
        if self.stop == 0:
            print("触发停机条件：达到最大迭代次数，注意此时可能未找到最优解，", end='')
            print("请增加迭代次数或增大梯度阈值重试。")
        elif self.stop == 1:
            print("触发停机条件：梯度小于阈值，运行正常")        
        
        return x

def _hesse(f, x):
    x_ = []
    for a in range(dim_max):
        x_.append(symbols('x'+str(a),real=True))  #设置符号变量
    f_ = f(x_)   # f的正则表达式f_
    dimensions = _dim(f_)
    
    item={}
    for i in range(dimensions):
        item[x_[i]] = x[i]

    dx = []     #求梯度 
    for a in range(dimensions):
        dx.append(diff(f_, symbols('x'+str(a),real=True))) #求偏导
    #print(f_)
    #print(dx)

    hesse = []  #求hesse阵
    for i in range(dimensions):
        ddx_i = []
        hesse_i = []
        for j in range(dimensions):
            ddx_i.append(diff(dx[i], symbols('x'+str(j),real=True)))
            #print(dx[i],'e',ddx_i[j])
            hesse_i.append(ddx_i[j].subs(item))
        hesse.append(hesse_i)
    return hesse