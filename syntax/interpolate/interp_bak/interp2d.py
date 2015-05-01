from libarray import *
from exceptions import Exception
from linalg import solve, det

def interp2d(x, y, data):

    class point(list):
        def __init__(self, *args):
	    list.__init__(self, args)
    pt0 = point(x, y)

    # check that the ref point is inside the data points
    xs, ys, fs = data.T
    if not(min(xs)<=x<=max(xs)) or not(min(ys)<=y<=max(ys)):
        raise Exception('Extrapolation is not allowed.')
    
    def dist(pt0, pt1):
        # compute the absolute distance
        return sqrt(sum([ (a-b)**2 for a, b in zip(pt0, pt1) ]))
    # generate a list of distance between the ref point and each data points
    dist0 = [ dist(point(x, y), pt0) for x, y, f in data ]
    
    if min(dist0)==0.0:
        indx = dist0.index(0.0)
	return data[indx][-1]

    A = array([])
    B = array([])
    for i in xrange(len(dist0)):
        # find the index of the nearest point
        indx = dist0.index(min(dist0))
	# replace the nearest point by infinite
        dist0[indx] = float('inf')
        x, y, f = data[indx]
	# build the system a*x+b*x+c = f -> linear interpolation
        A.append([x, y, 1.])
        B.append([f])
	# when the determinant is not null, the system can be solved
	det_A = det(A.T.dot(A))
	if det_A<>0:
	    break
    if det_A==0:
	raise Exception('The linear interpolation cannot be done.')

    # A.T.dot(A) gives a square matrix, independantly of the number of selected points 
    X = solve(A.T.dot(A), A.T.dot(B))
    a, b, c, = X.T[0]
    
    # check that the system give accurate result
    err = 0.0; ref = 0.0
    for [x, y, dummy], f in zip(A, B.T[0]):
        err += (a*x+b*y+c - f)**2
	ref += f*f
    err = sqrt(err/ref)
    if err>1e-6:
        Exception('The linear interpolation is not accurate.')
	return None
    
    return a*pt0[0]+b*pt0[1]+c

if __name__ == '__main__':
    x = array(range(-10,10,1), dtype=float)
    y = array(range(-10,10,1), dtype=float)

    data = array([])


    def f(x, y):
        return (x**2+y**2)
    
    for xi in x:
        for yi in y:
            data.append([xi, yi, f(xi, yi)])


    print interp2d(2.0001, 2, data)
     
