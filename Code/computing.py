
import time
from gurobipy import *
import json
import pandas as pd
import numpy as np


class data(object):
    def __init__(self):

        self.mark = 0
        self.mean_demand = {}
        self.lower_demand = None
        self.upper_demand = None
        self.random_demand = {}
        self.bar_x = None
        self.bar_z_j =None


    def initiData(self, *jsonFile_idx):

        jsonFile = r'/Users/fortunate/Downloads/pythonProject/EJOR-D-25-03831_R1/DataSet_9_21-5×4×3×3_20.json'

        with open(jsonFile) as f:
            res = json.load(f)
            f.close()

        '''baseParaDict input'''
        baseParaDict = res['baseParaDict']
        self.N = baseParaDict['N']
        self.T = baseParaDict['T']
        self.J = baseParaDict['J']
        self.K = baseParaDict['K']
        self.M = baseParaDict['M']
        self.L = baseParaDict['L']

        self.Tau=1
        self.zy = 0.8
        self.phi_k=1
        self.vartheta_p=1
        self.vartheta_r=1
        self.vartheta_h =1
        self.vartheta_f = 1



        cycle_len = 2
        self.TT = [i for i in range(self.T) if i % cycle_len == 0]
        self.TTT = [i for i in range(self.T) if i % cycle_len != 0]


        '''listParaDict Input'''
        listParaDict = res['listParaDict']
        self.c_jm = np.array(listParaDict['c_jm'])
        self.z_m = np.array(listParaDict['z_m'])
        self.b_l = np.array(listParaDict['b_l'])
        self.J_k = np.array(listParaDict['J_k'])
        self.K_j = np.array(listParaDict['K_j'], dtype=object)

        '''matrixParaDict Input'''
        matrixParaDict = res['matrixParaDict']
        self.l_jk = np.array(matrixParaDict['l_jk'])
        self.phi_jkm = np.array(matrixParaDict['phi_jkm'])
        self.f_jk = np.array(matrixParaDict['f_jk'])*self.vartheta_f
        self.h_nj = np.array(matrixParaDict['h_nj'])*self.vartheta_h
        self.r_nlj = np.array(matrixParaDict['r_nlj'])*self.vartheta_r
        self.p_n = np.array(listParaDict['p_n'])*self.vartheta_p

        def strIntotuple(str):
            '''change the str key into tuple'''
            str = str.strip('(')
            str = str.strip(')')
            splitStr = str.split(',')
            res = []
            for i in splitStr:
                res.append(int(i))
            res = tuple(res)
            return res


        '''Uniformaly distribution demand'''
        self.lower_demand = res['lower_demand']
        self.upper_demand = res['upper_demand']
        self.mean_demand=40

        for key, value in res['A_tilde_demand'].items():
            tuplekey = strIntotuple(key)
            data.random_demand[tuplekey] = value

data = data()
data.initiData(2) # 导入数据