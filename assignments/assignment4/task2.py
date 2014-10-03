# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 14:01:40 2014

@author: Axel
"""


from __future__ import division
from scipy import *
from matplotlib.pyplot import *

close("all")

def p(x):
    return x**9-18*x**8+144*x**7-672*x**6+2016*x**5-4032*x**4+5376*x**3-4608*x**2+2304*x-512

def p_(x):
    return (x-2)**9
    
def la(x):
    return x**9 - 9e1*x**8 + 36e2*x**7 - 84e3*x**6 + 126e4*x**5 - 126e5*x**4 + 84e6*x**3 - 36e7*x**2 + 9e8*x - 1e9
    
def la_(x):
    return (x-10)**9

xa = linspace(9.92, 10.08, 161)

plot(xa, la(xa))
figure()
plot(xa, la_(xa))

#x = linspace(1.92, 2.08, 161)
#
#plot(x,p(x))
#figure()
#plot(x, p_(x))


x = linspace(1.92, 2.08, 161)

plot(x,p(x))
figure()
plot(x,p_(x))
show('all')
