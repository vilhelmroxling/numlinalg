from __future__ import division
from scipy import *
from scipy import linalg as sla
from numpy import *
import numpy as np
from numpy import linalg as nla

def bisection(A):
    """
    A symmettric tridiagonal matrix.
    """
    n = A.shape[0]
    dets = zeros(n)
    counts = 0
    dets[0] = A[0, 0]
    if dets[0] < 0:
        counts += 1
    dets[1] = nla.det(A[:2,:2])
    if dets[0]*dets[1] < 0:
        counts += 1
    
    for i in range(2,n):
        dets[i] = A[i, i]*dets[i-1] - A[i, i-1]**2*dets[i-2]
        if dets[i]*dets[i-1] < 0:
            counts += 1

    return counts

def p(A, x):
    """
    A symmettric tridiagonal matrix.
    """
    n = A.shape[0]
    ps = zeros(n+2)
    counts = 0
    ps[-2] = 0
    ps[-1] = 1
    
    for i in range(n):
        ps[i] = (A[i, i]-x)*ps[i-1] - A[i, i-1]**2*ps[i-2]
        if ps[i]*ps[i-1] < 0:
            counts += 1

    return counts


def counteiginterval(A, a, b):
    n = A.shape[0]
    bcounts = p(A, b)
    acounts = p(A, a)
    return bcounts - acounts

def findeig(A, a, b, totcounts, eigs = [], tol = 1e-3):
    if b-a < tol:
        eigs.append([a, b])
    rightcounts = counteiginterval(A, (a+b)/2, b)
    leftcounts = totcounts - rightcounts
    if rightcounts > 0:
        findeig(A, (a+b)/2, b, rightcounts, eigs)
    if leftcounts > 0:
        findeig(A, a, (a+b)/2, leftcounts, eigs)

n = 10
a, b = 0, 1

A = random.rand(n, n)
A = (A + A.T)/2
A = sla.hessenberg(A)
eigs = []
findeig(A, -1, 1, counteiginterval(A, -1, 1), eigs)
print eigs
print sla.eig(A)[0]
