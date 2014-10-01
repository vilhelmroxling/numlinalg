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

x = linspace(1.92, 2.08, 161)

plot(x,p_(x))