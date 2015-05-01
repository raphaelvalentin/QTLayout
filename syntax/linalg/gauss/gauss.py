def gauss_pivot_total(A,B):
    # vecteur de pivotation des solutions

    a=[[j for j in i] for i in A]
    b=[j for j in B]

    n=len(b)
    pivot_sol=[0.]*n
    x=[0.]*n

    for i in xrange(0,n,1):
        pivot_sol[i]=i

    for k in xrange(0,n-1,1):
        # max pour le pivot total
        ref=0.0
        for i in xrange(k,n,1):
            for j in xrange(k,n,1):
                if ref<a[i][j]:
                    ref=a[i][j]
                    ligne=i
                    colonne=j
		elif ref<-a[i][j]:
                    ref=-a[i][j]
                    ligne=i
                    colonne=j

        # pivotations
        for j in xrange(k,n,1):
            a[k][j], a[ligne][j] = a[ligne][j], a[k][j]
        b[k], b[ligne] = b[ligne], b[k]
        
        for i in xrange(0,n,1):
            a[i][k], a[i][colonne] = a[i][colonne], a[i][k]

        # remplissage du vecteur accorde aux pivotations
        pivot_sol[k], pivot_sol[colonne] = pivot_sol[colonne], pivot_sol[k]

        if a[k][k]==0.0:
            return []

        # reduction
        for i in xrange(k+1,n,1):
            p=a[i][k]/a[k][k]
            for j in xrange(k,n,1):
                a[i][j] -= p*a[k][j]
            b[i] -= p*b[k]

    # resolution
    for i in xrange(n-1,-1,-1):
        s=0.0
        for j in xrange(i+1,n,1):
            s+=a[i][j]*b[j]
        b[i]=(b[i]-s)/a[i][i]

    # pivotation des solutions
    for i in xrange(0,n,1):
        x[pivot_sol[i]]=b[i]

    return x

if __name__ == '__main__':
    from random import random
    from time import time
    n = 200
    a = [[random() for i in xrange(n)] for i in xrange(n)]
    b = [[1.] for i in xrange(n)]
    t0 = time()
    x = wrap_gauss_pivot_total(a, array(b))
    from libarray import *
    x = array(x)
    A = array(a)
    
    print max(abs(array(A.dot(array(x)).T[0])-1.0))
