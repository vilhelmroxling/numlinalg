from    __future__  import division
from    scipy       import *
from    matplotlib.pyplot   import *

import numpy as np

class Orthogonalization:

    """
    It takes as an input the m x n array A (transposed for coding c
    onvenience) and assume_li which is a boolean(False by default) t
    elling us if the columns of A can be assumed to be linearly iden
    pendent or not. 
    """

    def __init__(self, A, assume_li = False):
        self.assume_li = assume_li
        self.A = A
        [self.n, self.m] = A.shape

    def gramschmidt(self):
        tol = 1e-12
        
        Q = np.zeros((self.n, self.m))
        Q[0] = self.A[0]/np.linalg.norm(self.A[0])
        lin_dep_inds = []

        for i in range(1,self.n):
            f = self.A[i]
            projvec = zeros(self.m)

            for k in range(i):
                projvec += np.dot(Q[k],f)*Q[k]        
        
            v = f - projvec
            vn = np.linalg.norm(v)

            if self.assume_li:
                Q[i] = v/vn

            else:
                if vn > tol:
                    Q[i] = v/vn
                else:
                    lin_dep_inds.append(i)
        
        return delete(Q, lin_dep_inds, axis=0)

    def householder(self, both = True):
        R = self.A.copy()
        self.v_list = []

        for i in range(self.n):
            a = R[i][i:]
            e = zeros(self.m-i)
            e[0] = np.linalg.norm(a)
            v = np.sign(a[0])*e + a
            v = v/np.linalg.norm(v)
            self.v_list.append(v)
            R[i,i:] = -sign(a[0])*e
            for j in range(1,self.n-i):
                R[i+j,i:] -= 2*v*np.dot(v, R[i+j,i:])
            
        if not both:
            return R, self.v_list
        
        Q = zeros((self.m, self.m))
        I = eye(self.m)
    
        for e, i in zip(I[:], range(self.m)):
            if self.n-1-i >= 0:
                k = self.n-1-i
                e[i:] -= 2*self.v_list[i][0]*self.v_list[i]
            else:
                k = 0
                e[self.n-1:] -= 2*self.v_list[-1][i-(self.n-1)]*self.v_list[-1]
                
            for j in range(k+1, self.n):
                v = self.v_list[self.n-1-j]
                e[self.n-1-j:] -= 2*np.dot(v, e[self.n-1-j:])*v

            Q[i] = e

        return Q.T, R.T
            
    def qdot(self, b):
        
        for i in range(self.n):
            v = self.v_list[self.n-1-i]
            b[self.n-1-i:] -= 2 * np.dot(v, b[self.n-1-i:])*v

        return b

    def givens(self, both = True):
        
        g_list = []
        R = self.A.copy()
        Q = eye(self.m)

        for i in range(self.n):
            
            if i == self.m-1:
                break

            #G = eye(self.m)
            x1 = R[i, -2]
            x2 = R[i, -1]
            r = np.sqrt(x1**2+x2**2)
            c = x1/r
            s = x2/r
            J = np.array([[c, s], [s, -c]])
            #G[-2:,-2:] = J
            R[:,-2:] = np.dot(R[:,-2:], J.T)
            if both:
                Q[-2:,:] = np.dot(J, Q[-2:,:])

            for j in range(3, self.m-i+1):
                x1 = R[i, -j]
                x2 = r
                r = sqrt(x1**2+x2**2)
                c = x1/r
                s = x2/r
                J = np.array([[c, s], [s, -c]])
                R[i:,-j:-j+2] = np.dot(R[i:,-j:-j+2], J.T)
                if both:
                    G[-j:-j+2,-j:] = np.dot(J,G[-j:-j+2,-j:])
                    #Q[-j:-j+2,:] = np.dot(J,Q[-j:-j+2,:])

            #if both:
             #   Q = np.dot(G, Q)

            g_list.append(G)
            #R = np.dot(R, G.T)

        if not both:
            #return g_list, R.T

        return Q.T, R.T
                
