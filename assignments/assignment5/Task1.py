# -*- coding: utf-8 -*-
"""
Created on Fri Oct 03 13:28:11 2014

@author: Axel
"""


from __future__ import division
from scipy import *
from matplotlib.pyplot import *



def switch(A, B):
    try:
        C, D = A.copy(), B.copy()
        return D,C
    except:
        return B, A

def invP(P):
    n = P.__len__()
    Pinv = []
    for k in range(n):
        i = P.index(k)
        Pinv.append(i)
    return Pinv
    
    
def LU(A):
    n = A.shape[0]
    U = A.copy()
    L = eye(n)
    p = range(n)
    
    for k in range(n-1):
        i = find(np.abs(U[k:,k]) == np.max(np.abs(U[k:,k]))) + k
        p[k], p[i] = switch(p[k],p[i])
        U[k,k:], U[i,k:] = switch(U[k,k:], U[i,k:])
        L[k,:k], L[i,:k] = switch(L[k,:k], L[i,:k])
        for j in range(k+1, n):
            ljk = U[j,k]/U[k,k]
            L[j,k] = ljk            
            U[j,k:] -= ljk*U[k,k:]
            
    return L, U, p
    
def fsub(L, b):
    n = b.__len__()
    x = zeros(n)
    
    for j in range (n):
        s = 0
        for i in range(j):
            s += L[j,i]*x[i]
        x[j] = (b[j] - s)/L[j,j]
        
    return x

def bsub(U, b):
    n = b.__len__()
    x = zeros(n)
    show = zeros((n,n))
    for j in range (n):
        s = 0
        for i in range(n-j,n):
            s += U[-j-1,i]*x[i]
            show[-j-1,i] = 1
            print show, "\n"
        x[-j-1] = (b[-j-1] - s)/L[-j-1,-j-1]
        
    return x


def solveMat(A,b):
    L, U, P = LU(A)
    pb = b[P]
    
    y = fsub(L,pb)
    x = bsub(U, y)
    
    return x
    

size = 5

A = random.rand(size,size)
b = random.rand(size)
L, U, P = LU(A)
x = solveMat(A,b)
print dot(A,x)-b

#print "______________________________", "\n"
#print P, "\n\n"
#print np.round(L, decimals=5), "\n\n"
#print np.round(U, decimals=5), "\n\n"
#print np.round(dot(L,U)-A[P], decimals = 5), "\n\n"
#
#print "______________________________", "\n\n"
#LU = dot(L,U)
#print np.round(LU[invP(P)] - A, decimals = 5)