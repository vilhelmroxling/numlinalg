# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 13:24:37 2014

@author: Axel
"""


from __future__ import division
from scipy import *
from scipy import linalg
from matplotlib.pyplot import *
import numpy as np

beta = 1e-5
alpha = 1
size = 50

#A = random.rand(size,size)
#Ainv = np.linalg.inv(A)

A = linalg.hilbert(size)
Ainv = linalg.invhilbert(size)

u,s,v = np.linalg.svd(A)

b = u[:,0]*alpha
db = u[:,-1]*beta

#x = np.linalg.solve(A,b)
#xh = np.linalg.solve(A,b + db)
x = np.dot(Ainv,b)
xh = np.dot(Ainv,b+db)
#dx = np.dot(Ainv, db)

dx = xh-x

k = np.linalg.norm(Ainv,ord=2)*np.linalg.norm(A,ord=2)
print k

vl = np.linalg.norm(dx, ord = 2)/np.linalg.norm(x, ord=2)

hl = (np.linalg.norm(db,ord=2)/np.linalg.norm(b,ord=2))
print vl/hl

