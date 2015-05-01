from libarray import *
from random import gauss
from function import *

def f(x):
    return x*x+2

def f(x):
    return (x-0.5)*(x-0.5)+2.0+1.*gauss(0, 0.00001)


def fmin(func, x0, x2, args=(), xtol=1e-4, ftol=1e-4, maxiter=100):
    fx0 = func(x0, *args)
    fx2 = func(x2, *args)
    mem = {'x':array([x0, x2]), 'f(x)':array([fx0, fx2]), 'solution':float('nan')}
    for i in xrange(maxiter):
        if fx2<fx0:
            x0, x2 = x2, x0
	    fx0, fx2 = fx2, fx0 
        x1 = 0.5*(x0+x2)
        fx1 = func(x1, *args)
	mem['x'].append(x1)
	mem['f(x)'].append(fx1)
        if fx1>0.0:
            x2 = x1
	    fx2 = fx1
        else:
            x0 = x1
	    fx0 = fx1
        if abs(fx1)<ftol or abs(x2-x0)<xtol:
            break
    v = closest1(mem['f(x)'], 0.0)
    i = mem['f(x)'].index(v)
    mem['solution'] = mem['x'][i]
    return mem

def fmin2(func, x0, x2, xtol=1e-4, ftol=1e-4, maxiter=50):
    fx0 = func(x0)
    fx2 = func(x2)
    mem = {'x':array([x0, x2]), 'f(x)':array([fx0, fx2]), 'solution':float('nan')}
    abserr = abs(mem['f(x)'])
    i = abserr.index(min(abserr))
    mem['solution'] = mem['x'][i]
    fx1 = float('nan')
    d = 1e300
    for i in xrange(maxiter):
        _fx1 = fx1
        x1 = 0.5*(x0+x2)
        fx1 = func(x1)
	mem['x'].append(x1)
	mem['f(x)'].append(fx1)
	v = closest1(mem['f(x)'], 0.0)
	i = mem['f(x)'].index(v)
	i0 = mem['f(x)'].index(max([i for i in mem['f(x)'] if i<0]))
    	x0 = mem['x'][i0]
	i2 = mem['f(x)'].index(min([i for i in mem['f(x)'] if i>0]))
	x2 = mem['x'][i2]
        if abs(fx1)<ftol or abs(x2-x0)<xtol:
            break
	try:
	    if abs(d-(mem['f(x)'][i0]-mem['f(x)'][i2]))<ftol:
	        break
	except:
	    pass    
	d = mem['f(x)'][i0]-mem['f(x)'][i2]    
    v = closest1(mem['f(x)'], 0.0)
    i = mem['f(x)'].index(v)
    mem['solution'] = mem['x'][i]
    return mem

	    
#r = fmin(f, -3, 5., maxiter=200, ftol=1e-6)
#if r:
#    print
#    print 'result', r['solution']
#    print 'iter', len(r['x'])
#    for i, v in enumerate(zip(r['x'], r['f(x)'])):
#        print i, v[0], v[1]


#from random import gauss
#def f(x):
#    x0 = x[0]
#    return abs(x0*x0-2.0+gauss(0, 0.1))


#from optimize.barecmaes2 import fmin
#from optimize.cma import fmin
#r = x = fmin(f, 1 * [10], 0.5, ftarget=1e-3)
#if r:
#    print
#    print 'result', r

#from math import sqrt
#x = [sqrt(2)]
#print [f(x) for i in xrange(10)]
