from math import *
from exceptions import Exception
from libarray import *

__all__ = ['svd', 'lstsq', 'SingularValueDecomposition']

_DEBUG = False

def shape(obj):
    _shape = []
    if isinstance(obj, list):
        _shape.append(len(obj))
        if isinstance(obj[0], list):
            _shape.extend(shape(obj[0]))
    return tuple(_shape)

def dot(mat1, mat2):
    m, n, p = len(mat1), len(mat2[0]), len(mat1[0])
    result = [[0.0 for j in xrange(n)] for i in xrange(m)]
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(p):
                result[i][j] += mat1[i][k]*mat2[k][j]
    return result
 
def transpose(mat):
    m, n = shape(mat)
    return [[mat[i][j] for i in xrange(m)] for j in xrange(n)]

def hypo(a, b):
    if abs(a) > abs(b):
        r = b/a
        r = abs(a)*sqrt(1.0 +  r*r)
    elif b != 0.0:
        r = a / b;
        r = abs(b) * sqrt(1.0 + r*r)
    else:
        r = 0.0
    return r
        

class SingularValueDecomposition:

    def __init__(self, Arg):
        """
        Construct the singular value decomposition 
        Derived from LINPACK code
        @param A Rectangular matrix
        @return Structure to access U, S and V
        """

        # Internal storage of U
        self.U = []
        # Internal storage of V
        self.V = []
        # Internal storage of singular values
        self.s = []
        # Row dimension
        self.m = 0
        # Column dimension
        self.n = 0

        # Initialize
        self.m, self.n = shape(Arg)
        A = []
        for i in xrange(self.m):
            row = []
            for j in xrange(self.n):
                row.append(Arg[i][j])
            A.append(row)
            
        self.s = [0. for i in xrange(self.n)]
        nu      = min(self.m, self.n)
        e       = [0.0 for i in xrange(self.n)]
        work    = [0.0 for i in xrange(self.m)]
        wantu   = True
        wantv   = True      
        nct = min(self.m - 1, self.n)
        nrt = max(0, min(self.n - 2, self.m))
        
        self.U = [[0.0 for i in xrange(min(self.m + 1, self.n))] for j in xrange(self.m)]
        self.V = [[0.0 for i in xrange(self.n)] for j in xrange(self.n)]

        # Reduce A to bidiagonal form, storing the diagonal elements
        # in s and the super-diagonal elements in e.
    
        for k in xrange(0, max(nct,nrt), 1):
    
            if k < nct:
                # Compute the transformation for the k-th column and
                # place the k-th diagonal in s[k].
                # Compute 2-norm of k-th column without under/overflow.
                self.s[k] = 0
                for i in xrange(k, self.m, 1):
                    self.s[k] = hypo(self.s[k], A[i][k])
                if (self.s[k] != 0.0): 
                    if A[k][k] < 0.0:
                        self.s[k] = -self.s[k]
                    for i in xrange(k, self.m, 1):
                        A[i][k] /= self.s[k]
                    A[k][k] += 1.0
                self.s[k] = -self.s[k]
                    
            for j in xrange(k + 1, self.n, 1):
                if (k < nct) and (self.s[k] != 0.0):
                    # Apply the transformation.
                    t = 0
                    for i in xrange(k, self.m, 1):
                        t += A[i][k] * A[i][j]
                    t = -t / A[k][k]
                    for i in xrange(k, self.m, 1):
                        A[i][j] += t * A[i][k]
                    # Place the k-th row of A into e for the
                    # subsequent calculation of the row transformation.
                    e[j] = A[k][j]
                    
    
            if (wantu and (k < nct)): 
                # Place the transformation in U for subsequent back
                # multiplication.
                for i in xrange(k, self.m, 1):
                    self.U[i][k] = A[i][k]
    
            if (k < nrt):
            # Compute the k-th row transformation and place the
            # k-th super-diagonal in e[k].
            # Compute 2-norm without under/overflow.
                e[k] = 0
                for i in xrange(k + 1, self.n, 1):
                    e[k] = hypo(e[k], e[i])
                if (e[k] != 0.0):
                    if (e[k+1] < 0.0):
                        e[k] = -e[k]
                    for i in xrange(k + 1, self.n, 1):
                        e[i] /= e[k]
                    e[k+1] += 1.0
            
                e[k] = -e[k]
                if ((k+1 < self.m) and (e[k] != 0.0)): 
                    # Apply the transformation.
                    for i in xrange(k+1, self.m, 1):
                        work[i] = 0.0
                    for j in xrange(k+1, self.n, 1):
                        for i in xrange(k+1, self.m, 1):
                            work[i] += e[j] * A[i][j]
                    for j in xrange(k + 1, self.n, 1):
                        t = -e[j] / e[k+1]
                        for i in xrange(k + 1, self.m, 1):
                            A[i][j] += t * work[i]
    
                if (wantv): 
                  # Place the transformation in V for subsequent
                  # back multiplication.
                  for i in xrange(k + 1, self.n, 1):
                      self.V[i][k] = e[i]
             
        # Set up the final bidiagonal matrix or order p.
        p = min(self.n, self.m + 1)
        if (nct < self.n):
            self.s[nct] = A[nct][nct]
        if (self.m < p):
            self.s[p-1] = 0.0
        if (nrt + 1 < p):
            e[nrt] = A[nrt][p-1]
        e[p-1] = 0.0
        
        # If required, generate U.
        if (wantu):
            for j in xrange(nct, nu, 1): 
                for i in xrange(0, self.m, 1):
                    self.U[i][j] = 0.0
                self.U[j][j] = 1.0
            for k in xrange(nct - 1, -1, -1):
                if (self.s[k] != 0.0): 
                    for j in xrange(k + 1, nu, 1): 
                        t = 0
                        for i in xrange(k, self.m, 1):
                            t += self.U[i][k] * self.U[i][j]
                        t = -t / self.U[k][k]
                        for i in xrange(k, self.m, 1):
                            self.U[i][j] += t * self.U[i][k]
                    for i in xrange(k, self.m, 1):
                        self.U[i][k] = -self.U[i][k]
                    self.U[k][k] = 1.0 + self.U[k][k]
                    for i in xrange(0, k - 1, 1):
                       self.U[i][k] = 0.0
                else: 
                    for i in xrange(0, self.m, 1):
                        self.U[i][k] = 0.0
                    self.U[k][k] = 1.0
        
        # If required, generate V.
        if (wantv):
            for k in xrange(self.n - 1, -1, -1):
                if (k < nrt) and (e[k] != 0.0):
                    for j in xrange(k + 1, nu, 1):
                        t = 0
                        for i in xrange(k + 1, self.n, 1):
                            t += self.V[i][k]* self.V[i][j]
                        t = -t / self.V[k+1][k]
                        for i in xrange(k + 1, self.n, 1):
                            self.V[i][j] += t * self.V[i][k]
              
                for i in xrange(0, self.n, 1):
                    self.V[i][k] = 0.0
                self.V[k][k] = 1.0
            
        # Main iteration loop for  xrange the singular values.
        pp   = p - 1
        iter = 0
        eps  = pow(2.0, -52.0)
        
        while p > 0:
              
            # Here is where a test for too many iterations would go.
            # This section of the program inspects for negligible
            # elements in the s and e arrays.  On completion the
            # variables kase and k are set as follows:
            # kase = 1  if s(p) and e[k-1] are negligible and k<p
            # kase = 2  if s(k) is negligible and k<p
            # kase = 3  if e[k-1] is negligible, k<p, and
            #           s(k), ..., s(p) are not negligible (qr step).
            # kase = 4  if e(p-1) is negligible (convergence).
    
            for k in xrange(p - 2, -2, -1): 
                if k == -1:
                    break
                if abs(e[k]) <= eps * (abs(self.s[k]) + abs(self.s[k+1])):
                    e[k] = 0.0
                    break
            
            if (k == p - 2):
                kase = 4
            else: 
                for ks in xrange(p - 1, k-1, -1): 
                    if ks == k:
                        break
                    t = (ks != p) * abs(e[ks]) + (ks != k + 1 ) * abs(e[ks-1])
                    if (abs(self.s[ks]) <= eps * t):  
                        self.s[ks] = 0.0
                        break
              
                if ks == k:
                    kase = 3
                elif (ks == p-1):
                    kase = 1
                else: 
                    kase = 2
                    k = ks
            k += 1 
            
            # Perform the task indicated by kase.
            if kase==1:
                # Deflate negligible s(p).
                f = e[p-2]
                e[p-2] = 0.0
                for j in xrange(p - 2, k-1, -1): 
                    t  = hypo(self.s[j],f)
                    cs = self.s[j] / t
                    sn = f / t
                    self.s[j] = t
                    if (j != k): 
                        f = -sn * e[j-1]
                        e[j-1] = cs * e[j-1]
                    if (wantv): 
                        for i in xrange(0, self.n, 1): 
                            t = cs * self.V[i][j] + sn * self.V[i][p-1]
                            self.V[i][p-1] = -sn * self.V[i][j] + cs * self.V[i][p-1]
                            self.V[i][j] = t
            
            # Split at negligible s(k).
            elif kase==2:
                f = e[k-1]
                e[k-1] = 0.0
                for j in xrange(k, p, -1): 
                    t = hypo(self.s[j], f)
                    cs = self.s[j] / t
                    sn = f / t
                    self.s[j] = t
                    f = -sn * e[j]
                    e[j] = cs * e[j]
                    if (wantu): 
                        for i in xrange(0, self.m, 1): 
                            t = cs * self.U[i][j] + sn * self.U[i][k-1]
                            self.U[i][k-1] = -sn * self.U[i][j] + cs * self.U[i][k-1]
                            self.U[i][j] = t
                   
            # Perform one qr step.
            elif kase==3:
                # Calculate the shift.                          
                scale = max(max(max(max(
                        abs(self.s[p-1]),abs(self.s[p-2])),abs(e[p-2])),
                        abs(self.s[k])), abs(e[k]))
                sp   = self.s[p-1] / scale
                spm1 = self.s[p-2] / scale
                epm1 = e[p-2] / scale
                sk   = self.s[k] / scale
                ek   = e[k] / scale
                b    = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0
                c    = (sp * epm1) * (sp * epm1)
                shift = 0.0
                if (b != 0.0) or (c != 0.0): 
                    shift = sqrt(b * b + c)
                    if b < 0.0:
                        shift = -shift
                    shift = c / (b + shift)
                f = (sk + sp) * (sk - sp) + shift
                g = sk * ek  
                # Chase zeros.  
                for j in xrange(k, p-1, 1): 
                    t  = hypo(f,g)
                    cs = f/t
                    sn = g/t
                    if j != k:
                        e[j-1] = t
                    f = cs * self.s[j] + sn * e[j]
                    e[j] = cs * e[j] - sn * self.s[j]
                    g = sn * self.s[j+1]
                    self.s[j+1] = cs * self.s[j+1]
                    if wantv: 
                        for i in xrange(0, self.n, 1): 
                            t = cs * self.V[i][j] + sn * self.V[i][j+1]
                            self.V[i][j+1] = -sn * self.V[i][j] + cs * self.V[i][j+1]
                            self.V[i][j] = t
     
                    t = hypo(f,g)
                    cs = f/t
                    sn = g/t
                    self.s[j] = t
                    f = cs * e[j] + sn * self.s[j+1]
                    self.s[j+1] = -sn * e[j] + cs * self.s[j+1]
                    g = sn * e[j+1]
                    e[j+1] = cs * e[j+1]
                    if wantu and (j < self.m - 1): 
                        for i in xrange(0, self.m, 1): 
                            t = cs * self.U[i][j] + sn * self.U[i][j+1]
                            self.U[i][j+1] = -sn * self.U[i][j] + cs * self.U[i][j+1]
                            self.U[i][j] = t
    
                e[p-2] = f
                iter = iter + 1
            
            # Convergence.
            elif kase==4:
                # Make the singular values positive.
                if self.s[k] <= 0.0: 
                    self.s[k] = (self.s[k] < 0.0) * -self.s[k]
                    if wantv: 
                        for i in xrange(0, pp+1, 1):
                            self.V[i][k] = -self.V[i][k]
                 
                # order the singular values.  
                while k < pp: 
                    if self.s[k] >= self.s[k+1]:
                        break
                    t = self.s[k]
                    self.s[k] = self.s[k+1]
                    self.s[k+1] = t
                    if wantv and (k < self.n - 1): 
                        for i in xrange(0, self.n, 1): 
                            t = self.V[i][k+1]
                            self.V[i][k+1] = self.V[i][k]
                            self.V[i][k] = t
                  
                
                    if wantu and (k < self.m-1): 
                        for i in xrange(0, self.m, 1): 
                            t = self.U[i][k+1]
                            self.U[i][k+1] = self.U[i][k]
                            self.U[i][k] = t
                  
                    k += 1
               
                iter = 0
                p -= 1
            # end switch      
        # end while        
        # end constructor
        
    def getU(self):
        """ 
        Return the left singular vectors
        @return U
        """
        #return self.U
        U = [[0.0 for i in xrange(min(self.m, self.n))] for j in xrange(self.m)]
        #U = [[0.0 for i in xrange(self.m)] for j in xrange(min(self.m+1, self.n))]
        for i in xrange(self.m):
            for j in xrange(min(self.m, self.n)):
                U[i][j] = self.U[i][j]
        return U

    def getV(self):
        """ 
        Return the right singular vectors
        @return V
        """
        return self.V
         
    def getSingularValues(self):
        """ 
        Return the one-dimensional array of singular values
        @return diagonal of S
        """
        return self.s

    def getS(self):
        """
        Return the diagonal matrix of singular values
        @return S
        """
        S = [[0.0 for i in xrange(self.n)] for j in xrange(self.n)]
        for i in xrange(0, self.n, 1):
            for j in xrange(0, self.n, 1):
                S[i][j] = 0.0
            S[i][i] = self.s[i]
        return S
        
    def norm2(self):
        """
        Two norm
        @return max(S)
        """
        return self.s[0]
        
    def cond(self):
        """
        Two norm condition number
        @return max(S)/min(S)
        """
        return self.s[0] / self.s[min(self.m, self.n) - 1]

    def rank(self):
        """
        Effective numerical matrix rank
        @return Number of nonnegligible singular values.
        """
        eps = pow(2.0, -52.0)
        tol = max(self.m, self.n) * self.s[0] * eps
        r = 0
        for i in xrange(0, len(self.s), 1):
            if (self.s[i] > tol):
              r+=1
        return r


def svd(A):
    if not(_DEBUG):
        m, n = shape(A)
	if m<n:
            raise Exception('Error on shape of the matrix: %d x %d Matrix is underestimated'%(m, n))
    UsV = SingularValueDecomposition(A)
    s = array(UsV.getSingularValues())
    U = array(UsV.getU())
    VT = array(transpose(UsV.getV()))
    return U, s, VT


def lstsq(A, B):
    # computing the inverse using the SVD decomposition
    U, s, V = svd(A)
    S_inv = [[0. for i in xrange(len(s))] for j in xrange(len(s))]
    for i in xrange(len(s)):
        if s[i]>pow(2, -52):
            S_inv[i][i] = 1./s[i]
    A_inv = dot(dot(transpose(V),S_inv),transpose(U))
    x = dot(A_inv,B) # solving Ax=b computing x = A^-1*b
    err = sqrt(sum((bi[0]-Axi[0])**2 for bi, Axi in zip(B, dot(A, x))))
    return x, err, sum(1 for x in s if x>0), s
    
def _TestCase1(n, m):
    from random import random
    import numpy as np
    a = [[random() for i in xrange(m)] for i in xrange(n)]
    b = [[1.] for i in xrange(n)]
    x1 = lstsq(a, b)[0]
    x2 = np.linalg.lstsq(a, b)[0]
    b1 = dot(a, x1)
    b2 = dot(a, x2)
    err1 = sqrt(sum([(bi[0][0]-bi[1][0])**2 for bi in zip(b1, b)]))
    err2 = sqrt(sum([(bi[0][0]-bi[1][0])**2 for bi in zip(b2, b)]))
    err = 2*(err1-err2)/(err1+err2)
    if abs(err)>3:
        raise Exception('Error on Test Case: %d x %d Matrix does not pass'%(n, m))



if __name__ == '__main__':
    _DEBUG = True
    try:
        print 'TestCase1 : square matrix 10x10...',
        for i in xrange(100):
            _TestCase1(10, 10)
        print 'pass'
    except:
        print 'does not pass'
    try:
        print 'TestCase2 : overdetermined matrix 15x5...',
        for i in xrange(100):
            _TestCase1(15, 5) 
        print 'pass'
    except:
        print 'does not pass'
    try:
        print 'TestCase3 : underdetermined matrix 5x15...',
        for i in xrange(100):
            _TestCase1(5, 15) 
        print 'pass'
    except:
        print 'does not pass'
    exit()


    """
    a = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    b = [[1237.], [1941.], [2417.]]
    U, s, V = svd(a)
    Ur, sr, Vr = np.linalg.svd(a)
    print U
    print Ur
    print s==sr
    print V==Vr
    exit()
    """
    
    a = transpose([[1,0,0],[0,1,0],[0,0,1],[-1,1,0],[-1,0,1],[0,-1,1]])
    b = [[1237.], [1941.], [2417.], [711.], [1177.], [475.]]
    
    a = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0., 0., 2.]]
    b = [[1237.], [1941.], [2417.], [4834.]]


    import numpy as np 
    U, s, V = np.linalg.svd(a)
    #print U
    #print s
    #print V.T
    #from libarray import array
    #print solve(a, b)
    print np.linalg.lstsq(a, b)
    #U, s, V = svd(a)
    print lstsq(a, b)
    exit(0)
    #print np.array(U)
    #print np.array(s)
    #print np.array(V).T
    S_inv = [[0. for i in xrange(len(s))] for j in xrange(len(s))]
    for i in xrange(len(s)):
        if s[i]>pow(2, -52):
            S_inv[i][i] = 1./s[i]
    #print transpose(V)
    a1 = dot(transpose(V),S_inv)
    a2 = transpose(U)
    #print np.array(a1)
    #print np.array(a2)
    A_inv = dot(dot(transpose(V),S_inv),transpose(U))
    x = dot(A_inv,b) # solving Ax=b computing x = A^-1*b
    print x

    exit()

    


    
    import numpy as np 
    U, s, V = np.linalg.svd(a)
    S_inv = [[0. for i in xrange(len(s))] for j in xrange(len(s))]
    for i in xrange(len(s)):
        if s[i]>pow(2, -51):
            S_inv[i][i] = 1./s[i]
   # computing the inverse using the SVD decomposition
    A_inv = dot(dot(V.T,S_inv),U.T)
    x = dot(A_inv,b) # solving Ax=b computing x = A^-1*b
    print x
    


