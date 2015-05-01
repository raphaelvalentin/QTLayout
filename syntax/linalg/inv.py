from libarray import *
from numpy import linalg

__all__ = ['inv']

def inv(X):
    Y = linalg.inv(X)
    return array([[x for x in row] for row in Y])
    
