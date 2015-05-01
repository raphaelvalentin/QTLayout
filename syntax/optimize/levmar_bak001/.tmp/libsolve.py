# -*- coding: cp1252 -*-
from math import exp, sqrt
from libarray import *

def norm(M):
    err = 0.
    for i in M:
        for j in i:
            err+=j*j
    return sqrt(err)

def zero(n):
    x=[]
    for i in xrange(n):
        x.append(0.)
    return x

def wrap_gauss_pivot_total(A,B):
    B = B.transpose()[0]
    X = gauss_pivot_total(A,B)
    Z=array([])
    for x in X:
        Z.append([x])
    return Z


def gauss_pivot_total(A,B):
    # vecteur de pivotation des solutions

    a=[[j for j in i] for i in A]
    b=[j for j in B]

    n=len(b)
    pivot_sol=zero(n)
    x=zero(n)

    for i in xrange(0,n,1):
        pivot_sol[i]=i

    for k in xrange(0,n-1,1):
        # max pour le pivot total
        ref=0.0
        for i in xrange(k,n,1):
            for j in xrange(k,n,1):
                if abs(a[i][j])>ref:
                    ref=abs(a[i][j])
                    ligne=i
                    colonne=j

        # pivotations
        for j in xrange(k,n,1):
            a[k][j], a[ligne][j] = a[ligne][j], a[k][j]
        b[k], b[ligne] = b[ligne], b[k]
        
        for i in xrange(0,n,1):
            a[i][k], a[i][colonne] = a[i][colonne], a[i][k]

        # remplissage du vecteur accorde aux pivotations
        pivot_sol[k], pivot_sol[colonne]= pivot_sol[colonne], pivot_sol[k]

        if a[k][k]==0:
            return []

        # réduction
        for i in xrange(k+1,n,1):
            p=a[i][k]/a[k][k]
            for j in xrange(k,n,1):
                a[i][j]=a[i][j]-p*a[k][j]
            b[i]=b[i]-p*b[k]

    # résolution
    for i in xrange(n-1,-1,-1):
        s=0.0
        for j in xrange(i+1,n,1):
            s=s+a[i][j]*b[j]
        b[i]=(b[i]-s)/a[i][i]

    # pivotation des solutions
    for i in xrange(0,n,1):
        x[pivot_sol[i]]=b[i]

    return x
