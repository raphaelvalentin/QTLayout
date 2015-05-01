from exceptions import Exception
from libarray import *

class LinAlgError(Exception):
    pass

import lu
import gauss
import qr
import svd
import inv
import chol

gauss = gauss
LULinAlgError = lu.LinAlgError
lu = lu.LUDecomposition
qr = qr.householder
lstsq = svd.lstsq
svd = svd.svd
inv = inv.inv
chol = chol.cholesky

__all__=['lu', 'solve', 'gauss', 'norm', 'det', 'linfit', 'qr', 'svd', 'lstsq', 'inv', 'chol']

def norm(vec):
    return sqrt(sum([x*x for x in flatten(vec)]))

def solve(A, B, verbose=False):
    """
    solve using LU or SVD decomposition method.
    """
    try:
        ShapeB = shape(B)
	if len(ShapeB)==1:
	    B = [[element] for element in B]
	elif not( len(ShapeB)==2 and ShapeB[1]==1 ):
            raise ValueError, 'Setting B with a correct shape'
    except:
        raise ValueError, 'Setting B element with a sequence'
    try:
        ShapeA = shape(A)
	if len(ShapeA)<>2:
	    raise ValueError
    except:
        raise ValueError, 'Setting A element with a sequence and a correct shape'
    if ShapeA[0]==ShapeA[1]:
        try:
	    if verbose: print 'Try using LU decomposition...',
	    LU = lu(A)
	    X = array(LU.solve(B))
	    if verbose: print 'Solved.'
	except LULinAlgError as details:
	    if verbose: 
	        print 'Not solved because %s.'%details
		print 'Try using SVD decomposition...',
	    X = lstsq(A, B)[0]
	    if verbose: print 'Solved.'
	except Exception as details:
	    if verbose: print 'Not solved because %s.'%details
	    raise Exception('Unexpected Error')
    else:
	if verbose: print 'The matrix is not square. Try using SVD decomposition...',
        X = lstsq(A, B)[0]
	if verbose: print 'Solved.',
    if len(ShapeB)==2 and ShapeB[1]==1:
        return X
    elif len(ShapeB)==1:
        return array([element[0] for element in X])

        
def det(M):
    """Compute the determinant of a square matrix by Gaussian elimination"""
    M = [ list(row) for row in M ]
    n = len(M)
    res = 1.0
    for j in xrange(n-1, 0, -1):
        pivot, i = max((abs(M[k][j]), k) for k in xrange(j+1))
        pivot = M[i][j]
        if pivot == 0.0:
            return 0.0
        M[i], M[j] = M[j], M[i]
        if i != j:
            res = -res
        res *= pivot
        fact = -1.0/pivot
        for i in xrange(j):
            f = fact * M[i][j]
            for k in xrange(j):
                M[i][k] += f * M[j][k]
    res *= M[0][0]
    return res

def linfit(data):
    #calcul de la regression lineaire de donnee en 2 colonnes xi,yi
    # retourne a, b et r
    N=len(data)
    # calcul de <x> et <y>
    mean_x=0.
    mean_y=0.
    for i in xrange(N):
        mean_x+=data[i][0]
        mean_y+=data[i][1]
    mean_x=mean_x/N
    mean_y=mean_y/N

    # calcul de a
    num_a=0.
    den_a=0.
    for i in xrange(N):
        num_a+=data[i][0]*data[i][1]
        den_a+=(data[i][0]-mean_x)**2
    num_a=num_a-N*mean_x*mean_y
    a=num_a/den_a
    # calcul de b
    b=mean_y-a*mean_x
    # calcul de r
    num_r=0.
    den_r=0.
    for i in xrange(N):
        num_r+=(data[i][0]-mean_x)**2
        den_r+=(data[i][1]-mean_y)**2
    r=a*sqrt(num_r/den_r)    

    return [a,b,r]

