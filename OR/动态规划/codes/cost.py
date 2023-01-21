import numpy as np
#from numpy import array

class solution:
    def __init__(self, rt, st, ct, max_age) -> None:
        '''
        rt:年收入,st:设备年龄age时残值,ct:年运营成本,pt:购入设备成本,max_age:设备最大年龄
        '''
        self.max_age = max_age  # 设备寿命
        self.rt = rt    # 年营业额
        self.st = st    # 残值
        self.ct = ct    # 年运营成本
        # self.pt = pt    # 每阶段购进新设备的价格
        if type(rt) is float:
            self.rt = [rt]*(max_age+1)
        if type(st) is float:
            self.st = [st]*(max_age+1)
        if type(ct) is float:
            self.ct = [ct]*(max_age+1)
        #if type(pt) is float:
        #    self.pt = [pt]*(max_age+1)
        pass

    def equipment_renew(self, pt, stop_stage, age=1):
        '''pt:购入设备成本, pro:最大利润, h:最优规划（用R表示更新，K表示继续使用，E表示结束）'''
        if type(pt) is float:
            self.pt = [pt]*(stop_stage+1)
        else:
            self.pt = pt

        self.stop_stage = stop_stage
        h, pro = self._update_stage(age, stage=1)
        h.reverse()
        return h, pro
        
    def _update_stage(self, age, stage):
        '''对于state阶段，年龄为age的设备的规划'''
        if(age > self.max_age):
            return ['err'], -1000.
        if(stage > self.stop_stage):
            return ['E'], self.st[age]

        h_K_next, pro_K_next = self._update_stage(age+1, stage+1)
        h_R_next, pro_R_next = self._update_stage(1, stage+1)
        pro_K = pro_K_next + self.rt[age] - self.ct[age]
        pro_R = pro_R_next + self.rt[0] - self.ct[0] + self.st[age] - self.pt[stage-1]
        h_K = h_K_next + ['K']
        h_R = h_R_next + ['R']

        if pro_K >= pro_R:
            return h_K, pro_K
        else:
            return h_R, pro_R

    def equipment_renew_advanced(self, pt, stop_stage, age=1):
        ''' 改进算法：递归 '''
        '''pt:购入设备成本, pro:最大利润, h:最优规划（用R表示更新，K表示继续使用，E表示结束）'''
        if type(pt) is float:   # 处理pt
            self.pt = [pt]*(stop_stage)
        else:
            self.pt = pt

        DP_Table = np.zeros([stop_stage+1, self.max_age+1])         # 记录不同阶段状态下的最优解
        DP_Table_step = np.zeros([stop_stage+1, self.max_age+1])    # 记录不同状态下的最优选择，K/R/E
        h = []                                                      # 记录最优解

        for t in range(self.max_age+1):
            DP_Table[stop_stage][t] = self.st[t]
            DP_Table_step[stop_stage][t] = -1                       # -1表示结束(E)

        for T in range(stop_stage-1, -1, -1):
            for t in range(1, self.max_age):
                pro_K = DP_Table[T+1][t+1] + self.rt[t] - self.ct[t]    # 选择保留
                pro_R = DP_Table[T+1][1] + self.rt[0] - self.ct[0] - self.pt[T] + self.st[t]    # 选择更新
                if pro_K > pro_R:
                    DP_Table[T][t] = pro_K
                    DP_Table_step[T][t] = 0                         # 0表示保留(K)
                else:
                    DP_Table[T][t] = pro_R
                    DP_Table_step[T][t] = 1                         # 1表示更新(R)

            DP_Table[T][self.max_age] = DP_Table[T+1][1] + self.rt[0] - self.ct[0] \
                - self.pt[T] + self.st[self.max_age]    # 设备寿命结束，更新
            DP_Table_step[T][self.max_age] = 1

        t = age
        pro = DP_Table[0][t]
        for T in range(stop_stage+1):
            tmp = DP_Table_step[T][t]
            if tmp == 0:
                h.append('K')
                t = t+1
            elif tmp == 1:
                h.append('R')
                t = 1
            else:
                h.append('E')
                break
        return h, pro
        #return DP_Table_step, DP_Table

    '''
    def equipment_renew_advanced_2(self, stop_stage, age=0):
        # 改进算法：转化为台阶问题，要求pt固定
        # 设备i年更换一次，使用期间的总收益计算
        pro_period = 

        # 转化为上台阶问题进行计算
        h, pro = self._update_stage_1(pro_period, stop_stage, age)
        h.reverse()
        return h, pro
    '''