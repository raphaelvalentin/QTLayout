from barecmaes2_wrap import *
from math import *
from libarray import *
from random import random
from function import linspace

xmeas = linspace(0, 1.2, step=0.1)
def f(x, u0, ua, ub, vth):
    tox = 1e-1
    return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
ymeas = array([f(x, 300, -0.04, 0.04, -0.20) for x in xmeas])

def objective_parallel(Xs):
    Y = []
    for X in Xs:
        y = f(xmeas, *X)
	Y.append( sqrt(sum( ((ymeas-y))**2 )) )
    return Y

def objective(X):
    y = f(xmeas, *X)
    return sqrt( sum( ((ymeas-y))**2 ) )


X = [500., 0.01, 0.01, 0.2]
opts = {
        'scaling_of_variables':[1e2, 1e-2, 1e-2, 1e-1], 
	'bounds':[(50, float('inf')), (-0.1, 0.1), (-0.1, 0.1), (-0.5, 0.5)],
	'tolfun': 1e-12,
	'tolfunhist':1e-14,
	'tolx': 1e-014, 
	'parallel': True,
	'verb_disp' : 100,
	#'seed':5, 
	#'popsize' : 20,
	}
res = fmin(objective_parallel, X, 1, **opts)
print 
print 'x_best:', res[1]
print 'f(x_best):', res[0]

