from libarray import *
from linalg import solve

class LinearNDInterpolator(object):
    def __init__(self, points, values):
        self.points = points
	self.values = values
        (self.n, self.m) = shape(self.points)
	self.ranges = [(min(column), max(column)) for column in points.T]
       
    def __call__(self, point):
        for c, (min_r, max_r) in zip(point, self.ranges):
	    if not(min_r<=c<=max_r):
	        return float('nan')
	
        # generate a list of normalized distance between the ref point and each data points
        distance = [ (sum((a-b)**2/(a*a+b*b) for a, b in zip(pt, point) if a<>b), indx) for indx, pt in enumerate(self.points) ]
	# sort  the distance
	distance.sort()
        
        dist, indx = distance[0]
        if dist==0.0:
	    return self.values[indx]
	
        # initialise
        # build the system a*x+b*y+c*z...+k = f -> linear interpolation
        A = array([])
        B = array([])

        # 1st point
        dist, indx = distance[0]
        A.append( list(self.points[indx]) + [1.0] )
        B.append([self.values[indx]])

        for i in xrange(self.m):
	    # ith point
            for dist, indx in distance[1:]:
	        if not(self.points[indx][i] in A.T[i]):
                    A.append( list(self.points[indx]) + [1.0] )
                    B.append([self.values[indx]])
	            break

        X = solve(A, B)
        r = 0.0
        for i in xrange(self.m):
            r += X[i][0]*point[i]
	r += X[-1][0]
	return r

     
if __name__ == '__main__':
    
    x = array(range(-15,15,1), dtype=float)/5.
    y = array(range(-15,15,1), dtype=float)/5.
    z = array(range(-15,15,1), dtype=float)/5.
    print len(x)*len(y)*len(z), 'points'

    points = array([])
    values = array([])


    def f(x, y, z):
        return (x**2+y**2+z**2)
    
    for xi in x:
        for yi in y:
            for zi in z:
	        points.append([xi, yi, zi,])
		values.append(f(xi, yi, zi))


    LI = LinearNDInterpolator(points, values)
    print LI([1.15, 1.52, 2.51])
    print f(*[1.15, 1.52, 2.51])
     
