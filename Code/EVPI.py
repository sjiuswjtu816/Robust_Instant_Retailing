# time: 2022-01-15 23:19
import datetime
import time
from computing import data
from gurobipy import *
import pandas as pd
import numpy as np


startTime = time.time()
solver_Time = []
onsalesrevenue_List = []
storesalesrevenue_List = []

fixedCost_List = []

replenishmentCost_List = []
fulfillmentCost_List = []
holdingCost_List = []
obj_List = []

class Deter_Model(object):

    def __init__(self):
        self.sample = 0
        self.model = Model('Deter Problem')
        self.solver_time = 0

    def updateSample(self, sample):
        self.sample = sample

    def creatVar(self):
        self.x = self.model.addVars(data.J,data.M ,vtype=GRB.BINARY, name='x_') # 是否建造商店
        self.y = self.model.addVars(data.J, data.K,data.M, vtype=GRB.BINARY, name='y_')   #这里y,zZ含义互换了，要注意
        self.zZ = self.model.addVars(data.K, data.M, vtype=GRB.BINARY, name='z_')


        self.g = self.model.addVars(data.J, data.K,data.M,lb=0, vtype=GRB.CONTINUOUS, name='g_') #MNL线性化辅助变量
        self.d = self.model.addVars(data.K, lb=0, vtype=GRB.CONTINUOUS, name='d_') #MNL线性化辅助变量
        self.dlt = self.model.addVars(data.N,data.T, data.L, data.J, vtype=GRB.BINARY, name='dlt_')
        self.v = self.model.addVars(data.N, data.T,data.L, data.J, lb=0, vtype=GRB.CONTINUOUS, name='v_') #补货数量
        self.w = self.model.addVars(data.N, data.T, data.J,data.K, lb=0, vtype=GRB.CONTINUOUS, name='w_') #线上履约量
        self.q = self.model.addVars(data.N, data.T, data.J, data.K, lb=0, vtype=GRB.CONTINUOUS, name='q_') #线下销售量
        self.z = self.model.addVars(data.N, data.T+1, data.J, lb=0, vtype=GRB.CONTINUOUS, name='z_') #仓库库存


    def setObjective(self):
        '''objective fuction'''
        onsalesrevenue = LinExpr(0)
        storesalesrevenue = LinExpr(0)
        fixedCost = LinExpr(0)

        replenishmentCost = LinExpr(0)
        fulfillmentCost = LinExpr(0)
        holdingCost = LinExpr(0)


        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        onsalesrevenue.addTerms(data.p_n[n],self.w[n,t,j,k])

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        storesalesrevenue.addTerms(data.p_n[n],self.q[n,t,j,k])

        for j in range(data.J):
            for m in range(data.M):
                fixedCost.addTerms(data.c_jm[j,m], self.x[j,m])




        for n in range(data.N):
            for t in range(data.T):
                for j in range(data.J):
                    for l in range(data.L):
                        replenishmentCost.addTerms(data.r_nlj[n,l, j], self.v[n, t,l,j])

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        fulfillmentCost.addTerms(data.f_jk[j,k],self.w[n,t,j,k])

        for n in range(data.N):
            for t in range(data.T):
                for j in range(data.J):
                    holdingCost.addTerms(data.h_nj[n, j], self.z[n, t+1,j])





        totalrevenue = onsalesrevenue+storesalesrevenue-fixedCost-replenishmentCost-fulfillmentCost-holdingCost

        self.model.setObjective(totalrevenue, sense=GRB.MAXIMIZE)

    def addConstr(self):

        for n in range(data.N):
            for j in range(data.J):
                self.model.addConstr(self.z[n,0,j]==0,name='c_13'+str([n,j]))

        for k in range(data.K):
            self.model.addConstr(quicksum(self.x[j,m] for j in data.J_k[k] for m in range(data.M))>=1, name='c_14' + str([k]))


        for j in range(data.J):
            self.model.addConstr(quicksum(self.x[j,m] for m in range(data.M))<=1, name='c_15' + str([ j]))

        for k in range(data.K):
            for m in range(data.M):
                self.model.addConstr(quicksum(self.x[j, m] for j in data.J_k[k]) <= 2*self.zZ[k,m],name='c_16' + str([k,m]))


        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.model.addConstr(self.y[j,k,m] <= self.x[j,m]  ,name='c_1' + str([j,k,m]))

        for k in range(data.K):
            for m in range(data.M):
                self.model.addConstr(quicksum(self.y[j,k,m] for j in data.J_k[k]) == self.zZ[k,m], name='c_2' + str([k,m]))

        for k in range(data.K):
            for k_ in range(data.K):
                for j in data.J_k[k]:
                    for m in range(data.M):
                        lhs = quicksum(self.y[j_,k,m] for j_ in data.J_k[k] if data.l_jk[j_,k] > data.l_jk[j,k])
                        lhs += quicksum(self.y[j_,k_,m] for j_ in data.J_k[k_] if (data.l_jk[j_,k] <= data.l_jk[j,k] and data.l_jk[j_,k_] > data.l_jk[j,k_]))
                        self.model.addConstr(lhs + quicksum(self.x[j,m] for m in range(data.M)) <= 1,name='c_3' + str([j,k,k_,m]))



        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    self.model.addConstr(self.z[n,t+1,j]==self.z[n,t,j]+quicksum(self.v[n,t,l,j] for l in range(data.L))-quicksum(self.q[n,t,j,k] for k in data.K_j[j] if data.K_j[j][0]>=0)-quicksum(self.w[n,t,j,k] for k in data.K_j[j] if data.K_j[j][0]>=0),name='c_6'+str([n,t,j]))


        for n in range(data.N):
            for t in data.TTT:
                for j in range(data.J):
                    self.model.addConstr(self.z[n,t+1,j]==self.z[n,t,j]-quicksum(self.q[n,t,j,k] for k in data.K_j[j] if data.K_j[j][0]>=0)-quicksum(self.w[n,t,j,k] for k in data.K_j[j] if data.K_j[j][0]>=0),name='c_6.'+str([n,t,j]))




        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        self.model.addConstr(self.q[n,t,j,k]<=quicksum(self.g[j,k,m]*data.phi_jkm[j,k,m]*data.random_demand[n, t, k][self.sample] for m in range(data.M)),name='c_7'+str([n,t,j,k]))



        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    mn=quicksum(self.g[j,k,m]*data.phi_jkm[j,k,m]*data.random_demand[n, t, k][self.sample] for m in range(data.M) for j in data.J_k[k])
                    self.model.addConstr(quicksum(self.w[n,t,j,k] for j in data.J_k[k])<=self.d[k]*data.random_demand[n,t,k][self.sample]+data.zy*mn-data.zy*quicksum(self.q[n,t,j,k] for j in data.J_k[k]),name='c_8'+str([n,t,k]))




        for k in range(data.K):
                self.model.addConstr(quicksum(data.phi_jkm[j_,k,m_]*self.g[j_,k,m_] for j_ in data.J_k[k] for m_ in range(data.M))+self.d[k]==1,name='c_9'+str([j,k]))

        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.model.addConstr(self.d[k] - self.g[j, k,m] <= 1 - self.y[j, k,m], name='c_10' + str([j, k,m]))

        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.model.addConstr(self.g[j, k,m] <= self.d[k], name='c_11' + str([j, k,m]))

        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.model.addConstr(self.g[j, k,m] <= self.y[j, k,m], name='c_12' + str([j, k,m]))



        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    self.model.addConstr(quicksum(self.dlt[n,t,l,j] for l in range(data.L)) -quicksum(self.x[j,m] for m in range(data.M))== 0, name='c_17' + str([n, j]))

        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    for l in range(data.L):
                        self.model.addConstr(self.v[n, t,l, j]  >= data.b_l[l]*self.dlt[n,t,l,j],name='c_18' + str([n,l, j]))

        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    for l in range(data.L):
                        self.model.addConstr(self.v[n, t,l, j]  <= data.b_l[l+1]*self.dlt[n,t,l,j],name='c_19' + str([n,l, j]))
        for t in data.TT:
             for j in range(data.J):
                self.model.addConstr(quicksum(self.z[n,t,j] for n in range(data.N))+quicksum(self.v[n, t,l, j] for n in range(data.N) for l in range(data.L) ) <= quicksum(data.z_m[m]*self.x[j,m] for m in range(data.M)), name='c_20' + str([ j]))



        self.model.setParam(GRB.Param.TimeLimit, 3600*4)
        #self.model.setParam('integralityFocus',1)
        self.model.setParam(GRB.Param.MIPGap, 0.005)
        self.model.setParam('Threads', 12)

    def chgCoef(self, sample):
        '''change the RHS of demand constraint'''
        self.sample = sample
        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        self.model.remove(self.model.getConstrByName(name='c_7'+str([n,t,j,k])))
                        self.model.addConstr(self.q[n,t,j,k]<=quicksum(self.g[j,k,m]*data.phi_jkm[j,k,m]*data.random_demand[n, t, k][self.sample] for m in range(data.M)),name='c_7'+str([n,t,j,k]))


        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    self.model.remove(self.model.getConstrByName(name='c_8' + str([n, t, k])))
                    mn = quicksum(self.g[j, k, m] * data.phi_jkm[j, k, m] * data.random_demand[n, t, k][self.sample] for m in range(data.M) for j in data.J_k[k])
                    self.model.addConstr(
                        quicksum(self.w[n, t, j, k] for j in data.J_k[k]) <= self.d[k] * data.random_demand[n, t, k][
                            self.sample] + data.zy * mn - data.zy * quicksum(self.q[n, t, j, k] for j in data.J_k[k]),
                        name='c_8' + str([n, t, k]))

    def solver(self, i):
        '''updated and optimzied'''
        start = time.time()
        if i >= 1:
            self.model.setParam('OutputFlag', 0)
        # self.model.setParam(GRB.Param.MIPGap, 0.01)
        self.model.update()
        self.model.optimize()
        self.solver_time = time.time()-start
        solver_Time.append(self.solver_time)

    def printAns(self, i):

        print("\ndeterministic Prob ObjVal: ", self.model.ObjVal, "solved Time: ", self.solver_time)

        actual_on_revenue = 0
        actual_store_revenue = 0
        actual_fixed_cost = 0

        actual_replenishment_cost = 0
        actual_fulfillment_cost = 0
        actual_holding_cost = 0

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        actual_on_revenue += data.p_n[n]*self.w[n,t,j,k].x

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        actual_store_revenue += data.p_n[n]*self.q[n,t,j,k].x

        for j in range(data.J):
            for m in range(data.M):
                actual_fixed_cost += data.c_jm[j,m]*self.x[j,m].x

        for n in range(data.N):
            for t in range(data.T):
                for l in range(data.L):
                    for j in range(data.J):
                        actual_replenishment_cost+=data.r_nlj[n,l,j]*self.v[n,t,l,j].x

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        actual_fulfillment_cost+=data.f_jk[j,k]*self.w[n,t,j,k].x


        for n in range(data.N):
            for t in range(data.T):
                for j in range(data.J):
                    actual_holding_cost+=data.h_nj[n,j]*self.z[n,t+1,j].x

        totalrevenue = actual_on_revenue + actual_store_revenue - actual_fixed_cost  - actual_replenishment_cost - actual_fulfillment_cost - actual_holding_cost
        fixedCost_List.append(actual_fixed_cost)
        replenishmentCost_List.append(actual_replenishment_cost)
        fulfillmentCost_List.append(actual_fulfillment_cost)
        holdingCost_List.append(actual_holding_cost)
        onsalesrevenue_List.append(actual_on_revenue)
        storesalesrevenue_List.append(actual_store_revenue)
        obj_List.append(totalrevenue)


        if self.sample == 0:

            print("0: fixed cost: ", round(actual_fixed_cost))
            print("1: replenishment cost: ", round(actual_replenishment_cost, 1))
            print("2: fulfillment cost: ", round(actual_fulfillment_cost, 1))
            print("3: holding cost: ", round(actual_holding_cost, 1))
            print("4: on revenue: ", round(actual_on_revenue, 1))
            print("5: store revenue: ", round(actual_store_revenue, 1))
            print("Total revenue: ", totalrevenue)


        def partI():
            usedFDC = []
            for j in range(data.J):
                if self.x[j].x >= 0.5:
                    usedFDC.append(j)


if __name__ == '__main__':
    startTime = time.time()
    a = ['UD']

    for co in range(1):
        if co == 0:
            print("-------------------------------Uniformly Distribution Demand--------------------------------------------")
        elif co == 1:
            solver_Time = []
            fixedCost_List = []
            replenishmentCost_List = []
            allocationCost_List = []
            FDCFulfillmentCost_List = []
            RDCFulfillmentCost_List = []
            RDCHoldingCost_List = []
            FDCHoldingCost_List = []
            DropShippingCost_List = []
            obj_List = []
            print("-------------------------------Normaly Distribution Demand--------------------------------------------")
            data.updateDemand_ND()
        else:
            solver_Time = []
            fixedCost_List = []
            replenishmentCost_List = []
            allocationCost_List = []
            FDCFulfillmentCost_List = []
            RDCFulfillmentCost_List = []
            RDCHoldingCost_List = []
            FDCHoldingCost_List = []
            DropShippingCost_List = []
            obj_List = []
            print("-------------------------------TwoPoint Distribution Demand--------------------------------------------")
            data.updateDemand_TP()
        DeterModel = Deter_Model()
        DeterModel.creatVar()
        DeterModel.setObjective()
        DeterModel.addConstr()
        start = time.time()
        DeterModel.solver(i=0)
        DeterModel.printAns(i=0)
        A = 500
        for i in range(1, A):
            print("=================sovleing sample {}=================".format(i))
            DeterModel.chgCoef(sample=i)
            DeterModel.solver(i=i)
            DeterModel.printAns(i=i)

        print("avg solved time : ", sum(solver_Time)/A)
        print("average Fixed Cost : ", sum(fixedCost_List)/A)
        print("average replenishment Cost : ", sum(replenishmentCost_List) / A)
        print("average fulfillment Cost : ", sum(fulfillmentCost_List)/A)
        print("average holding Cost : ", sum(holdingCost_List)/A)
        print("average on revenue : ", sum(onsalesrevenue_List)/A)
        print("average store revenue : ", sum(storesalesrevenue_List)/A)
        print("average Total revenue: ", sum(obj_List)/A)

        AllData = [solver_Time, fixedCost_List, replenishmentCost_List,  fulfillmentCost_List, holdingCost_List, onsalesrevenue_List, storesalesrevenue_List, obj_List]

        AllData = pd.DataFrame(AllData)

        fileName = r'EVPI-{}-{}×{}×{}×{}_sample-{}.xlsx'

        AllData.T.to_excel(fileName.format(a[co], data.N, data.T, data.J, data.K, A), header=["Solved_Time", "FixedCost","Replenishment",
                                             "fulfillment", "holding", "onrevenue", "storerevenue", "Totalrevenue"], index=False)