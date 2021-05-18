# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:04:25 2021

@author: Louis
"""
import matplotlib.pyplot as plt
import numpy as np
import math

col,row = 50,50

array1 = np.zeros((col,1))
array2 = np.zeros((col,1))
def cal_CovrianceMatrix(var):
    arrayA = np.zeros((col,row))
    for i in range(col):
        for j in range(row):
            arrayA[i][j]=(float)(var**(abs(i-j)))
    
    return arrayA

def cal_BhattacharyyaDistance(m1,m2,covMatrix1,covMatrix2):

    for i in range(col):
        array1[i] = m1
    for i in range(col):
        array2[i] = m2
    """
    np.transpose(x2_mean-x1_mean).dot( np.transpose((cov1-cov2)/2) ).dot( x2_mean-x1_mean )/8
    """
    meanDiff = ((1/8) * (array2- array1).T.dot(np.linalg.inv(0.5*(covMatrix1+covMatrix2)))).dot( array2- array1)
    print("mean difference :",meanDiff[0][0])
    temp  = np.linalg.det(covMatrix1)*np.linalg.det(covMatrix2)
    
    covDiff  = 0.5 * np.log(np.linalg.det((covMatrix1+covMatrix2)*0.5)/(temp**0.5))
    print("covariance difference: ",covDiff)
    ans = meanDiff[0][0]+covDiff
    print("BhattacharyyaDistance: ",ans)
    return ans
    
def cal_BhattacharyyaErrorBound(p1,p2,s,ans1):

    
    bayesError = ((p1*p2)**0.5)*math.exp(-ans1)
    
    
    print("Bayes error :" ,bayesError)
    
if __name__ == '__main__':
    m1,var1 = 0,0.9
    m2,var2 = 2.5,0.7
    p1,p2 = 0.5,0.5
    s = 0.5
    #CovMatrix
    covMatrix1 = cal_CovrianceMatrix(var1)
    covMatrix2 = cal_CovrianceMatrix(var2)
    
    
    
    #Carculate Bhawttacharyya Bound
    ans = cal_BhattacharyyaDistance(m1,m2,covMatrix1,covMatrix2)
    #Carculate Bhattacharyya error bound
    cal_BhattacharyyaErrorBound(p1,p2,s,ans)
 