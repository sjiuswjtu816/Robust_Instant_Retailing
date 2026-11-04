# time: 2022-02-28 20:13
import time
import sys
import pandas as pd
from computing import data
import json
import numpy as np
import matplotlib.pyplot as plt
from gurobipy import *

onsalesrevenue_List = []
storesalesrevenue_List = []
fixedCost_List = []
replenishmentCost_List = []
fulfillmentCost_List = []
holdingCost_List = []
obj_List = []
procurement_quantity_List = []
sales_volume_List = []
total_demand_List = []

Variable_x_List=[]
Variable_y_List=[]
Variable_z_List=[]

MPtime_list = []
SPtime_list = []

class SPM(object):

    def __init__(self):
        self._obj = 0
        self._starttime = time.time()
        self._solverTime = 0
        self.SP = Model('Subproblem')

    def SP_solver(self):

        self.x = self.SP.addVars(data.J, data.M, vtype=GRB.BINARY, name='x_')
        self.y = self.SP.addVars(data.J, data.K, data.M, vtype=GRB.BINARY, name='y_')
        self.zZ = self.SP.addVars(data.K, data.M, vtype=GRB.BINARY, name='z_')
        self.v = self.SP.addVars(data.N,data.T, data.L, data.J, lb=0, vtype=GRB.CONTINUOUS, name='v_')
        self.g = self.SP.addVars(data.J, data.K, data.M, lb=0, vtype=GRB.CONTINUOUS, name='g_')
        self.d = self.SP.addVars(data.K, lb=0, vtype=GRB.CONTINUOUS, name='d_')
        self.dlt = self.SP.addVars(data.N,data.T, data.L, data.J, vtype=GRB.BINARY, name='dlt_')
        self.A = self.SP.addVars(data.N, data.T, data.K, lb=0, vtype=GRB.CONTINUOUS, name='A_')
        self.G = self.SP.addVar(lb=0, vtype=GRB.CONTINUOUS, name='G_')
        self.E1 = self.SP.addVars(data.N, data.T, data.K, lb=0, vtype=GRB.CONTINUOUS, name='E1_')
        self.F1 = self.SP.addVars(data.N, data.T, data.K, lb=0, vtype=GRB.CONTINUOUS, name='F1_')

        for k in range(data.K):
            self.SP.addConstr(quicksum(self.x[j, m] for j in data.J_k[k] for m in range(data.M)) >= 1,
                                 name='c_14' + str([k]))

        for j in range(data.J):
            self.SP.addConstr(quicksum(self.x[j, m] for m in range(data.M)) <= 1, name='c_15' + str([j]))

        for k in range(data.K):
            for m in range(data.M):
                self.SP.addConstr(quicksum(self.x[j, m] for j in data.J_k[k]) <= 2 * self.zZ[k, m],
                                     name='c_16' + str([k, m]))


        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.SP.addConstr(self.y[j, k, m] <= self.x[j, m], name='c_1' + str([j, k, m]))

        for k in range(data.K):
            for m in range(data.M):
                self.SP.addConstr(quicksum(self.y[j, k, m] for j in data.J_k[k]) == self.zZ[k, m],
                                     name='c_2' + str([k, m]))

        for k in range(data.K):
            for k_ in range(data.K):
                for j in data.J_k[k]:
                    for m in range(data.M):
                        lhs = quicksum(self.y[j_, k, m] for j_ in data.J_k[k] if data.l_jk[j_, k] > data.l_jk[j, k])
                        lhs += quicksum(self.y[j_, k_, m] for j_ in data.J_k[k_] if
                                        (data.l_jk[j_, k] <= data.l_jk[j, k] and data.l_jk[j_, k_] > data.l_jk[j, k_]))
                        self.SP.addConstr(lhs + quicksum(self.x[j, m] for m in range(data.M)) <= 1,
                                             name='c_3' + str([j, k, k_, m]))

        for k in range(data.K):
            self.SP.addConstr(
                quicksum(data.phi_jkm[j_, k, m_] * self.g[j_, k, m_] for j_ in data.J_k[k] for m_ in range(data.M)) +
                self.d[k] == 1, name='c_9' + str([j, k]))

        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.SP.addConstr(self.d[k] - self.g[j, k, m] <= 1 - self.y[j, k, m],
                                         name='c_10' + str([j, k, m]))

        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.SP.addConstr(self.g[j, k, m] <= self.d[k], name='c_11' + str([j, k, m]))

        for k in range(data.K):
            for j in data.J_k[k]:
                for m in range(data.M):
                    self.SP.addConstr(self.g[j, k, m] <= self.y[j, k, m], name='c_12' + str([j, k, m]))

        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    self.SP.addConstr(quicksum(self.dlt[n,t, l, j] for l in range(data.L))-quicksum(self.x[j,m] for m in range(data.M))==0 , name='c_17' + str([n, j]))

        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    for l in range(data.L):
                        self.SP.addConstr(self.v[n, t,l, j] >= data.b_l[l] * self.dlt[n, t,l, j],
                                         name='c_18' + str([n, l, j]))

        for n in range(data.N):
            for t in data.TT:
                for j in range(data.J):
                    for l in range(data.L):
                        self.SP.addConstr(self.v[n, t,l, j] <= data.b_l[l + 1] * self.dlt[n, t,l, j],
                                         name='c_19' + str([n, l, j]))

        for n in range(data.N):
            for t in data.TTT:
                for j in range(data.J):
                    for l in range(data.L):
                        self.SP.addConstr(self.v[n, t,l, j]  == 0,
                                          name='c_199' + str([n, l, j]))

        for n in range(data.N):
            for t in data.TTT:
                for j in range(data.J):
                    for l in range(data.L):
                        self.SP.addConstr(self.dlt[n, t, l, j] == 0,
                                          name='c_191' + str([n, l, j]))


        Z = {} #Inventory level

        W = {} #Fulfillment quantity

        Q = {} #In-store sales quantity

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in range(data.J):
                        W[n, t, j, k] = {}
                        for subk in range(data.K):
                            for subt in range(t + 1):
                                W[n, t, j, k][subk, subt] = self.SP.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS,name='w_{} {}^{} {}[{}, {}]'.format(n, t, j,k, subk,subt))
                        W[n, t, j, k][0] = self.SP.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS,name='w_{} {}^{} {}, {}'.format(n, t, j, k, 0))

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in range(data.J):
                        Q[n, t, j, k] = {}
                        for subk in range(data.K):
                            for subt in range(t + 1):
                                Q[n, t, j, k][subk, subt] = self.SP.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS,name='q_{} {}^{} {}[{}, {}]'.format(n, t, j,k, subk,subt))
                        Q[n, t, j, k][0] = self.SP.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS,name='q_{} {}^{} {}, {}'.format(n, t, j, k, 0))

        for n in range(data.N):
            for t in range(data.T + 1):
                for j in range(data.J):
                    Z[n, t, j] = {}
                    if t > 0:
                        for k in range(data.K):
                            for subt in range(t):
                                Z[n, t, j][k, subt] = self.SP.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z_{}^{} {}, {} {}'.format(n, t, j, k, subt))
                    Z[n, t, j][0] = self.SP.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z_{}^{} {} {}'.format(n, t, j, 0))


        '''set objective function'''
        fixedCost = LinExpr(0)


        totalsalesrevenue = LinExpr(0)
        GG = self.G
        replenishmentCost = LinExpr(0)


        for j in range(data.J):
            for m in range(data.M):
                fixedCost.addTerms(data.c_jm[j,m], self.x[j, m])

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    totalsalesrevenue.addTerms(40,self.A[n,t,k])

        for n in range(data.N):
            for t in data.TT:
                for l in range(data.L):
                    for j in range(data.J):
                        replenishmentCost.addTerms(data.r_nlj[n,l,j],self.v[n,t,l,j] )

        totalrevenue = -fixedCost-replenishmentCost+totalsalesrevenue+GG

        self.SP.setObjective(totalrevenue, sense=GRB.MAXIMIZE)

        up=40
        lo=40

        c_1 = LinExpr(0)
        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        c_1.addTerms(data.p_n[n],W[n, t, j, k][0])
                        c_1.addTerms(data.p_n[n],Q[n, t, j, k][0])
                        c_1.addTerms(-data.f_jk[j,k], W[n, t, j, k][0])

        for n in range(data.N):
            for t in range(data.T):
                for j in range(data.J):
                    c_1.addTerms(-data.h_nj[n,j], Z[n, t+1, j][0])

        for n in range(data.N):
            for sigma in range(data.K):
                for tau in range(data.T):
                    c_1.addTerms(lo,self.E1[n,tau,sigma])
                    c_1.addTerms(-up, self.F1[n, tau, sigma])

        self.SP.addConstr(c_1 >= self.G,name='c1_RP_' )


        for n in range(data.N):
            for sigma in range(data.K):
                for tau in range(data.T):
                    lhh= quicksum(data.p_n[n]*W[n,t,j,k][sigma,tau]+data.p_n[n]*Q[n,t,j,k][sigma,tau]-data.f_jk[j,k]*W[n,t,j,k][sigma,tau] for k in range(data.K) for t in range(tau,data.T) for j in data.J_k[k])
                    ghh= quicksum(data.h_nj[n,j]*Z[n,t+1,j][sigma,tau] for j in range(data.J) for t in range(tau,data.T) )
                    self.SP.addConstr(lhh-ghh-self.E1[n,tau,sigma]+self.F1[n,tau,sigma]==self.A[n,tau,sigma],name='c1_CP_' + str([n, tau, sigma]))


        for n in range(data.N):
            for j in range(data.J):
                for t in data.TT:
                    self.SP.addConstr(Z[n, t + 1, j][0] == Z[n, t, j][0] +quicksum(self.v[n,t,l,j] for l in range(data.L)) - quicksum(Q[n,t,j,k][0] for k in data.K_j[j] if data.K_j[j][0]>=0) -  quicksum(W[n, t, j, k][0] for k in data.K_j[j] if data.K_j[j][0]>=0), name='c3_RP_' + str([n, t, j]))
                    for sigma in range(data.K):
                        for tau in range(t):
                            self.SP.addConstr(Z[n, t + 1, j][sigma, tau] == Z[n, t, j][sigma, tau]  - quicksum(Q[n,t,j,k][sigma,tau] for k in data.K_j[j] if data.K_j[j][0]>=0) - quicksum(W[n, t, j, k][sigma, tau] for k in data.K_j[j] if data.K_j[j][0]>=0),name='c3_CP_' + str([n, t, j]) + str([sigma, tau]))
                        self.SP.addConstr(Z[n, t + 1, j][sigma, t] +quicksum(Q[n,t,j,k][sigma,t] for k in data.K_j[j] if data.K_j[j][0]>=0) + quicksum(W[n, t, j, k][sigma, t] for k in data.K_j[j] if data.K_j[j][0]>=0) == 0,name='c3_CP_' + str([n, t, j]) + str([sigma, t]))


        for n in range(data.N):
            for j in range(data.J):
                for t in data.TTT:
                    self.SP.addConstr(Z[n, t + 1, j][0] == Z[n, t, j][0]  - quicksum(Q[n,t,j,k][0] for k in data.K_j[j] if data.K_j[j][0]>=0) -  quicksum(W[n, t, j, k][0] for k in data.K_j[j] if data.K_j[j][0]>=0), name='c3_RP_' + str([n, t, j]))
                    for sigma in range(data.K):
                        for tau in range(t):
                            self.SP.addConstr(Z[n, t + 1, j][sigma, tau] == Z[n, t, j][sigma, tau]  - quicksum(Q[n,t,j,k][sigma,tau] for k in data.K_j[j] if data.K_j[j][0]>=0) - quicksum(W[n, t, j, k][sigma, tau] for k in data.K_j[j] if data.K_j[j][0]>=0),name='c3_CP_' + str([n, t, j]) + str([sigma, tau]))
                        self.SP.addConstr(Z[n, t + 1, j][sigma, t] +quicksum(Q[n,t,j,k][sigma,t] for k in data.K_j[j] if data.K_j[j][0]>=0) + quicksum(W[n, t, j, k][sigma, t] for k in data.K_j[j] if data.K_j[j][0]>=0) == 0,name='c3_CP_' + str([n, t, j]) + str([sigma, t]))


        E2 = {}
        F2 = {}

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        E2[n, t, j,k] = self.SP.addVars(data.K, t+1, lb=0, vtype=GRB.CONTINUOUS,
                                                      name='E1_{} {} {}'.format(n, t, j))
                        F2[n, t, j,k] = self.SP.addVars(data.K, t+1, lb=0, vtype=GRB.CONTINUOUS,
                                                      name='F1_{} {} {}'.format(n, t, j))

                        self.SP.addConstr(Q[n, t, j, k][0]-quicksum(E2[n, t, j,k][sigma, tau] * lo for sigma in range(data.K) for tau in range(t+1)) +quicksum(F2[n, t, j,k][sigma, tau] * up for sigma in range(data.K) for tau in range(t+1))
                                          <= 0, name='c4_RP_'+str([n, t, j, k]))
                        for tau in range(t + 1):
                            for sigma in range(data.K):
                                if tau == t and sigma == k:
                                    self.SP.addConstr(Q[n,t,j,k][sigma,tau]+E2[n, t, j,k][sigma, tau] - F2[n, t, j,k][sigma, tau] == quicksum(self.g[j,k,m]*data.phi_jkm[j,k,m] for m in range(data.M)) ,name='demandConstr_' + str([n, t,j, k]))
                                else:
                                    self.SP.addConstr(Q[n,t,j,k][sigma,tau]+E2[n, t, j,k][sigma, tau] - F2[n, t, j,k][sigma, tau] == 0)


        E3 = {}
        F3 = {}

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    E3[n, t,k] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS,
                                                  name='E1_{} {} {}'.format(n, t, j))
                    F3[n, t, k] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS,
                                                  name='F1_{} {} {}'.format(n, t, j))

                    self.SP.addConstr(data.zy*quicksum(Q[n, t, j, k][0] for j in data.J_k[k])+quicksum(W[n,t,j,k][0] for j in data.J_k[k])-quicksum(E3[n, t,k][sigma, tau] * lo for sigma in range(data.K) for tau in range(t+1)) +quicksum(F3[n, t,k][sigma, tau] * up for sigma in range(data.K) for tau in range(t+1))<=  0,
                                          name='c5_RP_' + str([n, t, k]))
                    for tau in range(t + 1):
                        for sigma in range(data.K):
                            if tau == t and sigma == k:
                                self.SP.addConstr(data.zy*quicksum(Q[n, t, j, k][sigma, tau] for j in data.J_k[k]) +quicksum(W[n,t,j,k][sigma,tau] for j in data.J_k[k])+E3[n, t,k][sigma, tau] - F3[n, t,k][sigma, tau]   == self.d[k]+data.zy*quicksum(data.phi_jkm[j,k,m]*self.g[j,k,m] for j in data.J_k[k] for m in range(data.M)),name='demandConstr_' + str([n, t, k]))
                            else:
                                self.SP.addConstr(
                                        data.zy*quicksum(Q[n, t, j, k][sigma, tau] for j in data.J_k[k]) +quicksum(W[n,t,j,k][sigma,tau] for j in data.J_k[k])+E3[n, t,k][sigma, tau] - F3[n, t,k][sigma, tau]   ==  0)


        E8 = {}
        F8 = {}

        for j in range(data.J):
            for t in data.TT:
                E8[t , j] = self.SP.addVars(data.N,data.K, t , lb=0, vtype=GRB.CONTINUOUS,
                                                      name='E5_{}^{} {}'.format(n, t + 1, j))
                F8[t , j] = self.SP.addVars(data.N,data.K, t , lb=0, vtype=GRB.CONTINUOUS,
                                                      name='F5_{}^{} {}'.format(n, t + 1, j))

                self.SP.addConstr(quicksum(Z[n,t,j][0] for n in range(data.N))+quicksum(self.v[n,t,l,j] for n in range(data.N) for l in range(data.L))-quicksum( data.z_m[m] * self.x[j, m] for m in range(data.M))-quicksum(
                    E8[t , j][cps,sigma, tau] * lo for cps in range(data.N) for sigma in range(data.K) for tau in range(t)) + quicksum(
                    F8[t , j][cps,sigma, tau] * up for cps in range(data.N) for sigma in range(data.K) for tau in range(t)) <= 0, name='c77_RP_' + str([t, j]))
                self.SP.addConstrs((Z[cps, t , j][sigma, tau]  +E8[t , j][cps,sigma, tau] -F8[t , j][cps,
                        sigma, tau]  == 0
                                        for cps in range(data.N) for sigma in range(data.K) for tau in range(t)),
                                       name='c77_CP_' + str([t, j]))



        E5 = {}
        F5 = {}

        for n in range(data.N):
            for j in range(data.J):
                self.SP.addConstr(Z[n, 0, j][0] == 0, name='c7_RP_0')
                for t in range(data.T):
                    E5[n, t + 1, j] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS, name='E5_{}^{} {}'.format(n, t+1, j))
                    F5[n, t + 1, j] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS, name='F5_{}^{} {}'.format(n, t+1, j))


                    self.SP.addConstr(Z[n, t + 1, j][0] +quicksum(E5[n, t+1,j][sigma, tau] * lo for sigma in range(data.K) for tau in range(t+1)) -quicksum(F5[n, t+1,j][sigma, tau] * up for sigma in range(data.K) for tau in range(t+1))>=  0, name='c7_RP_'+str([n, t+1]))
                    self.SP.addConstrs((Z[n, t + 1, j][sigma, tau] -E5[n, t + 1, j][sigma, tau] + F5[n, t + 1, j][sigma, tau] == 0
                                        for sigma in range(data.K) for tau in range(t + 1)), name='c7_CP_'+str([n, t+1]))

        E6 = {}
        F6 = {}

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        E6[n, t, j, k] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS, name='E6_{} {}^{} {}'.format(n, t, j, k))
                        F6[n, t, j, k] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS, name='F6_{} {}^{} {}'.format(n, t, j, k))


                        self.SP.addConstr(W[n, t, j, k][0] +quicksum(E6[n, t,j,k][sigma, tau] * lo for sigma in range(data.K) for tau in range(t+1)) -quicksum(F6[n, t,j,k][sigma, tau] * up for sigma in range(data.K) for tau in range(t+1))>=  0, name='c8_RP_'+str([n, t, j, k]))
                        self.SP.addConstrs((W[n, t , j,k][sigma, tau] - E6[n, t , j,k][sigma, tau] + F6[n, t , j,k][
                            sigma, tau]  == 0
                                            for sigma in range(data.K) for tau in range(t + 1)),
                                           name='c8_CP_' + str([n, t ,j,k]))


        E7 = {}
        F7 = {}

        for n in range(data.N):
            for t in range(data.T):
                for k in range(data.K):
                    for j in data.J_k[k]:
                        E7[n, t, j, k] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS,
                                                         name='E7_{} {}^{} {}'.format(n, t, j, k))
                        F7[n, t, j, k] = self.SP.addVars(data.K, t + 1, lb=0, vtype=GRB.CONTINUOUS,
                                                         name='F7_{} {}^{} {}'.format(n, t, j, k))

                        self.SP.addConstr(Q[n, t, j, k][0] + quicksum(
                            E7[n, t, j, k][sigma, tau] * lo for sigma in range(data.K) for tau in
                            range(t + 1)) - quicksum(
                            F7[n, t, j, k][sigma, tau] * up for sigma in range(data.K) for tau in
                            range(t + 1))  >= 0, name='c9_RP_' + str([n, t, j, k]))
                        self.SP.addConstrs((Q[n, t, j, k][sigma, tau] - E7[n, t, j, k][sigma, tau] + F7[n, t, j, k][
                            sigma, tau]  == 0
                                            for sigma in range(data.K) for tau in range(t + 1)),
                                           name='c9_CP_' + str([n, t, j, k]))





        self.SP.setParam('OutputFlag', 1)
        self.SP.setParam('Threads', 12)
        self.SP.setParam('TimeLimit', 14400)
        self.SP.setParam('MIPGap', 0.005)
        #self.SP.setParam('LogFile', 'gurobi_log.txt')
        #self.SP.setParam('Method',2)
        #self.SP.setParam('NodeMethod', 2)


        # self.SP.Params.Crossover = 0
        # self.SP.Params.Aggregate = 0
        # self.SP.Params.BarConvTol = 1e-8

        self.SP.optimize()

        self._obj = self.SP.ObjVal
        self._solverTime = time.time() - self._starttime

        def calVar(sample,p,v):


            w = np.zeros(shape=(data.N, data.T, data.J, data.K), dtype=float)
            q = np.zeros(shape=(data.N, data.T, data.J, data.K), dtype=float)


            for n in range(data.N):
                for t in range(data.T):
                    for k in range(data.K):
                        for j in range(data.J):
                            w[n, t, j, k] = W[n, t, j, k][0].x + sum(W[n, t, j, k][key].x * data.random_demand[n, key[1], key[0]][sample] for key in W[n, t, j, k].keys() if key != 0)#+sum(W1[n, t, j, k][key].x * data.randomtheta1_demand[n, key[1], key[0]][sample] for key in W1[n, t, j, k].keys() if key != 0)+sum(W2[n, t, j, k][key].x * data.randomtheta2_demand[n, key[1], key[0]][sample] for key in W2[n, t, j, k].keys() if key != 0)
                            if w[n,t,j,k] <0:
                                w[n,t,j,k]=0

            for n in range(data.N):
                for t in range(data.T):
                    for k in range(data.K):
                        for j in range(data.J):
                            q[n, t, j, k] = Q[n, t, j, k][0].x + sum(
                                Q[n, t, j, k][key].x * data.random_demand[n, key[1], key[0]][sample] for key in
                                Q[n, t, j, k].keys() if key != 0)#+sum(Q1[n, t, j, k][key].x * data.randomtheta1_demand[n, key[1], key[0]][sample] for key in Q1[n, t, j, k].keys() if key != 0)+sum(Q2[n, t, j, k][key].x * data.randomtheta2_demand[n, key[1], key[0]][sample] for key in Q2[n, t, j, k].keys() if key != 0)
                            if q[n,t,j,k] < 0:
                                q[n,t,j,k]=0

            return w,q


        def calExpVar(p):

            actual_ = 0
            actual = self.G.x

            actual_fixed_cost = 0

            actual_replenishment_cost = 0

            for n in range(data.N):
                for t in range(data.T):
                    for k in range(data.K):
                        actual_+=40*self.A[n,t,k].x

            for j in range(data.J):
                for m in range(data.L):
                    actual_fixed_cost += data.c_jm[j,m] * p[j,m]

            for n in range(data.N):
                for t in range(data.T):
                    for l in range(data.L):
                        for j in range(data.J):
                            actual_replenishment_cost += data.r_nlj[n,l, j] * self.v[n,t,l,j].x

            totalrevenue =  - actual_fixed_cost - actual_replenishment_cost +actual_+actual

            print('Exp Total revenue: ', totalrevenue)


        if ans == False:
            print("SPM1 exp obj : ", SPM1._obj)
            p=np.zeros(shape=(data.J,data.M), dtype=int)

            for j in range(data.J):
                for m in range(data.M):
                    p[j,m]=self.x[j,m].x

            calExpVar(p)

            gam = np.zeros(shape=(data.J,data.K,data.M), dtype=float)
            for k in range(data.K):
                for j in data.J_k[k]:
                    for m in range(data.M):
                        gam[j,k,m] = self.g[j,k,m].x

            v1=np.zeros(shape=(data.N,data.T,data.L,data.J), dtype=float)
            for n in range(data.N):
                for t in range(data.T):
                    for j in range(data.J):
                        for l in range(data.L):
                            v1[n,t,l,j]=self.v[n,t,l,j].x
                            # print(f'v{[n, t, l, j]}=', v1[n, t, l, j])


            dd=np.zeros(shape=(data.K), dtype=float)
            for k in range(data.K):
                dd[k]=self.d[k].x

            y1 = np.zeros(shape=(data.J, data.K,data.M), dtype=int)
            for k in range(data.K):
                for j in data.J_k[k]:
                    for m in range(data.M):
                        y1[j,k,m]=self.y[j, k, m].x

            zZ1 = np.zeros(shape=( data.K, data.M), dtype=int)
            for k in range(data.K):
                for m in range(data.M):
                    zZ1[k,m]=self.zZ[ k, m].x


            A = 500
            for i in range(A):
                w,q = calVar(sample=i,p=p,v=v1)


    def Dispose(self):
        self.SP.dispose()

    def print(self, lb, ub):
        print("SubProblem:          Optimal Objective:{}         Solved Time:{}".format(self._obj, round(self._solverTime, 1)))
        print("\nNew Bounds：({}, {})".format(lb, ub))




if __name__ == '__main__':

    print("begin time: ", time.ctime())
    with open('Log_RO_{}-{}-{}-{}.txt'.format(data.N,data.T,data.J,data.K), 'w') as f_txt:
        print("begin time: ", time.ctime())
        f_txt.write("begin time : {}".format(time.ctime()))

        ans = False


        SPM1 = SPM()
        SPM1.SP_solver()
        print("terminate time: ", time.ctime())
        f_txt.write("Terminate time : {}".format(time.ctime()))
        f_txt.close()


