from exceptions import Exception

class LinAlgError(Exception):
    pass

class LUDecomposition(object):
    _LU = list()
    _m = 0
    _n = 0
    _pivsign = 0
    _piv = list()
    #
    #    Decompose from A a LU matrix
    #
    def __init__(self, A):
        if isinstance(A, list):
            # Use a "left-looking", dot-product, Crout/Doolittle algorithm.
            self._m = len(A)
            self._n = len(A[0])
            self._LU = [[A[i][j] for j in xrange(self._n)] for i in xrange(self._m)]
            self._piv = [i for i in xrange(self._m)]
            self._pivsign = 1
            LUrowi = LUcolj = [0. for i in xrange(self._m)]
            # Outer loop.
            for j in xrange(self._n):
                # Make a copy of the j-th column to localize references.
                for i in xrange(self._m):
                    LUcolj[i] = self._LU[i][j]
                # Apply previous transformations.
                for i in xrange(self._m):
                    LUrowi = self._LU[i]
                    # Most of the time is spent in the following dot product.
                    kmax = min(i, j)
                    s = 0.0
                    for k in xrange(kmax):
                        s += LUrowi[k] * LUcolj[k]
                    LUcolj[i] -= s
                    LUrowi[j] = LUcolj[i]
                # Find pivot and exchange if necessary.
                p = j
                for i in xrange(j+1, self._m, 1):
                    if abs(LUcolj[i]) > abs(LUcolj[p]):
                        p = i
                if (p != j):
                    for k in xrange(self._n):
                        t = self._LU[p][k]
                        self._LU[p][k] = self._LU[j][k]
                        self._LU[j][k] = t
                    k = self._piv[p]
                    self._piv[p] = self._piv[j]
                    self._piv[j] = k
                    self._pivsign = self._pivsign * -1
                # Compute multipliers.
                if (j < self._m) and (self._LU[j][j] != 0.0):
                    for i in xrange(j+1, self._m, 1):
                        self._LU[i][j] /= self._LU[j][j]
        else:
            raise TypeError, 'A is not an instance of list'
            
    #
    #    Extract from LU the L matrix
    #
    def L(self):
        L = [[0. for i in xrange(self._n)] for i in xrange(self._m)]
        for i in xrange(self._m):
            for j in xrange(self._n):
                 if i > j:
                     L[i][j] = self._LU[i][j]
                 elif i == j:
                     L[i][j] = 1.0
        return L
    L = property(fget=L, doc='Extract from LU the L matrix')      

    #
    #    Extract from LU the U matrix
    #
    def U(self):
        U = [[0. for i in xrange(self._n)] for i in xrange(self._n)]
        for i in range(self._n):
            for j in xrange(self._n):
                if i <= j:
                    U[i][j] = self._LU[i][j]
        return U
    U = property(fget=U, doc='Extract from LU the U matrix')        

    #
    #     Is the matrix nonsingular?
    #
    def isNonsingular(self):
        for j in xrange(self._n):
            if self._LU[j][j] == 0.:
                return False
        return True

    #
    #     Count determinants
    # 
    def det(self):
        if self._m == self._n:
            d = self._pivsign
            for j in xrange(self._n):
                d *= self._LU[j][j]
            return d
        else:
            raise LinAlgError('MatrixDimensionException')

    #
    #    Solve A*X = B
    #
    def solve(self, B):
        if len(B) == self._m:
            if self.isNonsingular():
                # Copy right hand side with pivoting
                nx = len(B[0])
                X  = [[B[self._piv[i]][j] for j in xrange(nx)] for i in xrange(self._m)]
                # Solve L*Y = B(piv,:)
                for k in xrange(self._n):
                    for i in xrange(k+1, self._n, 1):
                        for j in xrange(nx):
                            X[i][j] -= X[k][j] * self._LU[i][k]
                # Solve U*X = Y;
                for k in xrange(self._n-1, -1, -1):
                    for j in xrange(nx):
                        X[k][j] /= self._LU[k][k]
                    for i in xrange(k):
                        for j in xrange(nx):
                            X[i][j] -= X[k][j] * self._LU[i][k]
                return X
            else:
                raise LinAlgError('MatrixSingularException')
        else:
            raise LinAlgError('MatrixSquareException')

if __name__ == '__main__':
    from random import random
    from time import time
    lu = LUDecomposition
    A = [[1.2, 2.3],[-3.4, 4.5],[5., -6.]]
    B = [[1.5],[2.],[-0.5]]
    x = lu(A)
    print 'L', x.L
    print 'U', x.U
    print 'X', x.solve(B)


    n = 30
    a = [[random() for i in xrange(n)] for i in xrange(n)]
    b = [[1.] for i in xrange(n)]
    t0 = time()
    LU = lu(a)
    x = LU.solve(b)
    from libarray import *
    x = array(x)
    A = array(a)
    
    print max(abs(array(A.dot(x).T[0])-1.0))
    
    
    print 'time', time()-t0
            
        
