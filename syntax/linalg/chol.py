from math import sqrt
from libarray import *

class LinAlgError(Exception):
    pass


class cholesky(array):
    def __init__(self, A):
        n = len(A)
        array.__init__(self, [[0.0 for j in xrange(n)] for i in xrange(n)])
        try:
            for i in xrange(n):
                for k in xrange(i+1):
                    s = sum(self[i][j] * self[k][j] for j in xrange(k))
	    	    self[i][k] = sqrt(A[i][i] - s) if (i==k) else \
		              (1.0 / self[k][k] * (A[i][k]  - s))
        except ValueError:
            raise LinAlgError('Matrix is not positive definite - Cholesky decomposition cannot be computed')
    def solve(self, b):
        n = len(b)
        c = [0.0 for i in xrange(n)]
        for i in xrange(n):
            c[i] = (b[i] - sum( self[i][k] * c[k] for k in xrange(i) )) / self[i][i]
        x = [0.0 for i in xrange(n)]
        for i in xrange(n-1, -1, -1):
            x[i] = (c[i] - sum( self[k][i] * x[k] for k in xrange(i+1, n) )) / self[i][i]
        return x
        	   
        
if __name__ == '__main__':
    
    A = array([[2, 1, 1, 3, 2],
               [1, 2, 2, 1, 1],    
               [1, 2, 9, 1, 5],    
               [3, 1, 1, 7, 1],    
               [2, 1, 5, 1, 8]], dtype=float)
    b = array([1, 1, 1, 1, 1], dtype=float)
    L = cholesky( A )

    X =  L.solve(b)
    from linalg import norm
    from libarray import *

    print 'Must be near 0: ', norm(A.dot(array([X]).T)- b)

