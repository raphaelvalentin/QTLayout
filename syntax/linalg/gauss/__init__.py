from libarray import *
from gauss import gauss_pivot_total

__all__ = ['gauss', 'gauss_pivot_total']

def gauss(A, B):
    try:
        ShapeB = shape(B)
	if len(ShapeB)==2 and ShapeB[1]==1:
	    B = [element[0] for element in B]
	elif len(ShapeB)==1:
	    pass
	else:
            raise ValueError, 'Setting B with a correct shape'
    except:
        raise ValueError, 'Setting B element with a sequence'
    X =  gauss_pivot_total(A, B)
    if len(ShapeB)==2 and ShapeB[1]==1:
        return array([[element] for element in X])
    return X
        
