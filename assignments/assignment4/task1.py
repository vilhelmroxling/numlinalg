# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 12:51:53 2014

@author: Axel
"""


from __future__ import division
from scipy import *
from matplotlib.pyplot import *

close("all")
mlist = []
mlist.append(2)

#Full rank matrices
for i in range(1,8):
    mlist.append(mlist[i-1]*2)

if(False):
    for m in mlist:
        figure()
        for n in range(0,500):
            A = np.random.randn(m,m)/np.sqrt(m)
            lamb = np.linalg.eigvals(A)
            slamb = sort(lamb)
            nA= np.linalg.norm(A,ord=2)
            u,s,v = np.linalg.svd(A)
            

            plot(s[-1], '*')

            
print np.linalg.norm(A, ord=2)            
#Triangular matrices      
if(True):
    for m in mlist:
        figure()
        for n in range(0,200):
            A = np.random.randn(m,m)/np.sqrt(m)
            A = np.triu(A)
            lamb = np.linalg.eigvals(A)
            slamb = sort(lamb)
            nA= np.linalg.norm(A,ord=2)
            u,s,v = np.linalg.svd(A)
            
            plot(s[-1], '*')

show('all')
