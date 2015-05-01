from math import cos, pi, acos
from numpy import *
from random import random
from optimize.cma import cma

noise = array([(random()*2-1)*1e-3 for i in xrange(31)])
xmeas = array([-0.5 ,-0.45 ,-0.4 ,-0.35 ,-0.3 ,-0.25 ,-0.2 ,-0.15 ,-0.1 ,-0.05 ,0.0 ,0.05 ,0.1 ,
                0.15 ,0.2 ,0.25 ,0.3 ,0.35 ,0.4 ,0.45 ,0.5 ,0.55 ,0.6 ,0.65 ,0.7 ,0.75 ,0.8 ,0.85 ,0.9 ,0.95 ,1.0])
def f(x, u0, ua, ub, vth):
    tox = 1e-1
    return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
ymeas = array([f(x, 300, 0.04, 0.04, 0.20) for x in xmeas])
ymeas += noise



def boundary(a, b):
    a = float(a)
    b = float(b)
    def f(x): 
        return a + 0.5*(b-a) * ( 1.0 - cos( x*pi/(b-a) ) )
    def finv(y): 
        return (b-a)/pi * acos( 1.0 - 2.0*(y-a)/(b-a) )
    return f, finv

def boundaries(bounds, scaling_of_variables=None):
    _bounds = []
    if scaling_of_variables:
        for a, b, s in zip(bounds[0], bounds[1],scaling_of_variables):
            _bounds.append( boundary(a/s, b/s) )
    else:
        for a, b in zip(bounds[0], bounds[1]):
            _bounds.append( boundary(a, b) )
    return _bounds

def f(X):
    [u0, ua, ub, vth] = X
    tox = 1e-1 
    y = u0/(1.0+(ua/tox)*(xmeas-vth)+(ub/tox/tox)*(xmeas-vth)**2)
    return sum( ((ymeas-y))**2 )

X = [500., 0.01, 0.01, 0.2]
opts = {
        #'seed':1234, 
	'scaling_of_variables':[1e2, 1e-2, 1e-2, 1e-1], 
	'bounds':[[50, 0, 0, 0],[600, 0.1, 0.1, 0.5]],
	'tolfun': 1e-12,
	'tolfunhist':1e-12,
	'tolx': 1e-06, 
	'parallel': False,
	}

def fmin(f, x0, sigma=1, **kwargs):

    scaling_of_variables = kwargs.pop('scaling_of_variables', None)
    bounds = kwargs.pop('bounds', None)
    parallel = kwargs.pop('parallel', True)
    
    if bounds:
        bounds = boundaries(bounds, scaling_of_variables)
    
    def wrapper(func, scaling_of_variables=None, bounds=None):
        if bounds:
	   if scaling_of_variables:
	       def decorator(X):
	           Y = [ s*f(x) for x, (f, finv), s in zip(X, bounds, scaling_of_variables) ]
		   return func(Y)
	   else:
	       def decorator(X):
	           Y = [ f(x) for x, (f, finv) in zip(X, bounds) ]
		   return func(Y)
        else:
	   if scaling_of_variables:
	       def decorator(X):
	           Y = [ s*x for x, s in zip(X, scaling_of_variables) ]
		   return func(Y)
	   else:
	       decorator = func
        return decorator
    
    if scaling_of_variables:
	x0 = [ x/s for x, s in zip(x0, scaling_of_variables) ]
    if bounds:
	x0 = [ _finv(x) for x, (_f, _finv) in zip(x0, bounds) ]
    if scaling_of_variables or bounds:
        f = wrapper(f, scaling_of_variables, bounds)

    opts = cma.CMAOptions()
    opts.update( kwargs )
    es = cma.CMAEvolutionStrategy(x0, 1, opts)
    while not es.stop(): 
        X = es.ask()
	if parallel:
	    Y = f(X)
	else:
	    Y = [ f(x) for x in X]
        es.tell(X, Y)
        es.disp()
	
    if bounds:
        for i, (f, finv) in enumerate(bounds):
            es.best.__dict__['x'][i] = f(es.best.__dict__['x'][i])
    if scaling_of_variables:
        for i, s in enumerate(scaling_of_variables):
            es.best.__dict__['x'][i] *= s
        
    print 'Termination by ', 
    for k, v in  es.stop().iteritems():
        print k, '=', v,
    print
    #cma.pprint(es.best.__dict__)
    return es.best.__dict__['f'], list(es.best.__dict__['x'])
    
res = fmin(f, X, 1, **opts)
print 'res', res
#print 'f(x)', f(x)



