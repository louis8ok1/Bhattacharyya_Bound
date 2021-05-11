# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:04:25 2021

@author: Louis
"""
import matplotlib.pyplot as plt
import numpy as np


col,row = 50,50
arrayA = np.zeros((col,row))
array1 = np.zeros((col,1))
array2 = np.zeros((col,1))
def cal_CovrianceMatrix(var):
    
    for i in range(col):
        for j in range(row):
            arrayA[i][j]=(float)(var**(abs(i-j)))
    
    return arrayA

def cal_BhattacharyyaBound(m1,m2,covMatrix1,covMatrix2,p1,p2,s):
    for i in range(col):
        array1[i] = m1
    for i in range(col):
        array2[i] = m2
    
   
    
    ans = (((s)**2)/2) * (array2- array1).T.dot(np.linalg.inv(s*covMatrix1 + s*covMatrix1)).dot( array2- array1)
    temp = 0.5*np.log(np.linalg.det((s*covMatrix1)+(s*covMatrix2))/(np.linalg.det(covMatrix1)**s * np.linalg.det(covMatrix2))**s) 
    ans +=temp
    
    #檢查是否e^(-ans)在0.4~0.6間
    BhattacharyyaBound = np.exp(-ans)
    print(BhattacharyyaBound[0][0])
    #第二小題了
    #upperbound = np.sqrt(p1*p2) * np.exp(ans)
    #print(upperbound[0])
    return BhattacharyyaBound[0][0]
    
    

if __name__ == '__main__':
    m1,var1 = 0,0.9
    m2,var2 = 2.5,0.7
    p1,p2 = 0.5,0.5
    s1,s2 = 0.5,0.66
    #CovMatrix
    covMatrix1 = cal_CovrianceMatrix(var1)
    covMatrix2 = cal_CovrianceMatrix(var2)
    #print(covMatrix1.shape)
    
    
    #Carculate Bhawttacharyya Bound
    BhattacharyyaBound1 = cal_BhattacharyyaBound(m1,m2,covMatrix1,covMatrix2,p1,p2,s1)
    BhattacharyyaBound2 = cal_BhattacharyyaBound(m1,m2,covMatrix1,covMatrix2,p1,p2,s2)
    
    #Carculate Upper Bound
    upperbound1 = np.sqrt(p1*p2) * np.exp(BhattacharyyaBound1)
    print(upperbound1)
    upperbound2 = np.sqrt(p1*p2) * np.exp(BhattacharyyaBound2)
    print(upperbound2)
    