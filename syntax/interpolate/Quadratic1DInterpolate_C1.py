from libarray import *
from linalg import solve
import numpy as np

class Quadratic1DInterpolate:
    # Continuous-1
    def __init__(self, x, y):
        self.x = x
        self.y = y
   
        self.xscale = 1.0/(x[-1]-x[0]), -1.0/(x[-1]-x[0])*x[0]
        self.yscale = 1.0/(y[-1]-y[0]), -1.0/(y[-1]-y[0])*y[0]
    
        x = [xi*self.xscale[0] + self.xscale[1] for xi in x]
        y = [yi*self.yscale[0] + self.yscale[1] for yi in y]
        dim = (len(self.x)-1)*3
        a = zeros((dim,dim))
        b = zeros((dim,1))
        
        j = 0
        for i in xrange(len(x)-1):
            x0, y0 = x[i], y[i]
            a[j][i*3] = x0**2
            a[j][i*3+1] = x0
            a[j][i*3+2] = 1.0
            j = j+1
            
            x1, y1 = x[i+1], y[i+1]
            a[j][i*3] = x1**2
            a[j][i*3+1] = x1
            a[j][i*3+2] = 1.0
            j = j+1
            
        for i in xrange(1, len(x)-1, 1):
            x1, y1 = x[i], y[i]
            a[j][(i-1)*3] = 2.0*x1
            a[j][(i-1)*3+1] = 1.0
            a[j][i*3] = -2.0*x1
            a[j][i*3+1] = -1.0
            j = j+1
            
        a[j][0] = 1.0
        a[j][3] = -1.0
        
        j=0
        for i in xrange(len(x)-1):
            b[j][0] = y[i]
            j = j+1
            b[j][0] = y[i+1]
            j = j+1
        
        self.p = [x[0] for x in solve(a, b)]
         
    def __call__(self, x):
        for i in xrange(len(self.x)-1):
            if self.x[i]<=x<=self.x[i+1]:
                break
        x = x*self.xscale[0] + self.xscale[1]
        y = x**2*self.p[3*i+0] + x*self.p[3*i+1] + self.p[3*i+2]
        return (y-self.yscale[1])/self.yscale[0]




if __name__ == '__main__':

    x = array([-28., -10.0, 0.0, 18.0, 25.0])
    y = array([-32.0, -7.0, -5.0, 2.0, 3.0])

    x = array([-30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0 , 15.0 , 20.0, 25.0])
    y = array([-32.7334924416, -27.7404558739,-22.7532568334, -17.7772363956, -12.8231836489,  -7.91469765844, -3.10922686697, 1.39928758096, 5.25687094185, 8.85522316541, 11.9612031633, 14.6927575271])

    x = array([0.9999, 0.997, 0.99, 0.98, 0.97, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0])
    x.reverse()
    y = array([-1e-6, -0.07740155, -0.14106736, -0.198997487, -0.243104916, -0.3122499, -0.435889894, -0.526782688, -0.6, -0.714142843, -0.8, -0.866025404, -0.916515139, -0.953939201, -0.979795897, -1.2])
    y.reverse()


    f = Quadratic1DInterpolate(x, y)
    from function import *
    from plotting import plot

    plt = plot()
    plt.scatter(x, y)
    
    x = linspace(min(x), max(x), step=0.001)
    y = array([f(xi) for xi in x])
    #x, y = derivative(x, y)
    plt.plot(x, y)

    plt.savefig('myplt.png', force=True)

