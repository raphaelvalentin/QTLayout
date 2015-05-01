from libarray import *
from linalg import solve


class Quadratic1DInterpolate(object):
    def __init__(self, x, y):
        self.x = list(x)
	self.y = list(y)
    def __call__(self, x):
        if not(min(self.x)<=x<=max(self.x)):
	    return float('nan')
	
	distance = [((self.x[i]-x)**2, i) for i in xrange(len(self.x))]
	distance.sort()
        
        dist, indx = distance[0]
        if dist==0.0:
	    return self.y[indx]
	
        A = []
        B = []
	
        # 1st point
        dist, indx = distance[0]
        A.append( [self.x[indx]**i for i in xrange(2, -1, -1)] )
        B.append( [self.y[indx]] )
	
        # 2st point
        dist, indx = distance[1]
        A.append( [self.x[indx]**i for i in xrange(2, -1, -1)] )
        B.append( [self.y[indx]] )
	
        # 3st point
        dist, indx = distance[2]
        A.append( [self.x[indx]**i for i in xrange(2, -1, -1)] )
        B.append( [self.y[indx]] )
	
        X = solve(A, B)
	
	return sum(X[i][0]*x**(2-i) for i in xrange(3))

if __name__ == '__main__':

    from function import *
    from plotting import plot

    x = array([-28., -10.0, 0.0, 18.0, 25.0, 35])
    y = array([-32.0, -7.0, -5.0, 2.0, 3.0, 7.0])
    f = Quadratic1DInterpolate(x, y)
    plt = plot()
    plt.scatter(x, y)
    x = linspace(min(x), max(x), step=0.1)
    y = array([f(xi) for xi in x])
    #x, y = derivative(x, y)
    plt.plot(x, y)

    plt.savefig('myplt.png', force=True)
