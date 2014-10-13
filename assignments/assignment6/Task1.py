from __future__ import division
from scipy import *
from scipy import linalg as sla
from numpy import *
import numpy as np
from numpy import linalg as nla
import pdb
import time
import math

def tridiag(Ah):
    """
    Transforms the real symmetric, square matrix A to Hessenberg form, i.e it tridiagonalizes A.
    """
    A = Ah.copy()
    n = A.shape[0]
    for k in range(n-2):
        x = A[k+1:, k]
        e1 = zeros(n-(k+1))
        e1[0] = 1
        v = sign(x[0])*nla.norm(x)*e1 + x
        v = v/nla.norm(v)
        A[k+1:,k:] = A[k+1:,k:] - 2*outer(v, dot(v, A[k+1:,k:]))
        A[k:,k+1:] = A[k:,k+1:] - 2*outer(dot(A[k:,k+1:], v), v)
    return A


def QReig(Ah, count, tri = True, zerotol = 1e-15):
    """
    A real symmetric, square matrix. Returns a matrix of the same size with the eigenvalues on the diagonal. count is a variable that counts the number of QR iterations.
    """
    A = Ah.copy()
    n = A.shape[0]
    if n == 1:
        return A, count
    if tri:
        A = tridiag(A)
    while 1:
        count = count + 1
        mu = A[-1,-1]
        #mu = A[math.trunc(n/2), math.trunc(n/2)]
        Q, R = sla.qr(A - mu*eye(n))
        A = dot(R,Q) + mu*eye(n)
        sdiag = A[range(1,n), range(n-1)]
        if min(abs(sdiag)) < zerotol:
            indv = np.where(abs(sdiag) < zerotol)[0]
            for i, ind in zip(range(indv.__len__()), indv):
                A[ind+1, ind] = A[ind, ind+1] = 0
                if i == 0:
                    ind0 = 0
                else:
                    ind0 = indv[i-1]+1
                A[ind0:ind+1,ind0:ind+1], count = QReig(A[ind0:ind+1, ind0:ind+1], count)
            A[indv[-1]+1:, indv[-1]+1:], count = QReig(A[indv[-1]+1:, indv[-1]+1:], count)
            # ind = np.where(abs(sdiag) == min(abs(sdiag)))[0][0]
            # A[ind+1, ind] = A[ind, ind+1] = 0
            # A[:ind+1,:ind+1], count = QReig(A[:ind+1, :ind+1], count)
            # A[ind+1:, ind+1:], count = QReig(A[ind+1:, ind+1:], count)
            break
    return A, count

n = 5

A = 10*random.rand(n, n)
A = (A+A.T)/2
# v = random.rand(n)
# v = v/nla.norm(v)
# A = eye(n) - 2*outer(v, v)
t0 = time.time()
D, counts = QReig(A, 0, zerotol = 1e-6)
t1 = time.time()-t0
eigs = sla.eig(A)[0]
#print D, "\n", np.round(D, decimals = 3), "\n", eigs, "\n", nla.norm(sort(eigs)-sort(diag(D))), "\n", counts
print nla.norm(sort(eigs)-sort(diag(D))), "\n", counts, "\ntime: ", t1 

# t0 = time.time()
# D, counts = QReig(A, 0, tri = False)
# t1 = time.time()-t0
# eigs = sla.eig(A)[0]
# print nla.norm(sort(eigs)-sort(diag(D))), "\n", counts, "\ntime: ", t1 
