import barecmaes2 as cma
from barecmaes2 import *
from math import cos, pi, acos, sqrt, isnan, isinf
import random
from optimize.transformation import BoxContraintsTransformation, ScalingTransformation


__all__ = ['fmin']

class CMAES(cma.CMAES):
    def __init__(self, *args, **kwargs):
        self.tolx = kwargs.pop('tolx', 1e-11)
	self.tolfunhist = kwargs.pop('tolfunhist', 1e-12)
	kwargs['ftarget']  = kwargs.pop('tolfun', 1e-12)
        cma.CMAES.__init__(self, *args, **kwargs)
    def stop(self):
        """return satisfied termination conditions in a dictionary like 
        {'termination reason':value, ...}, for example {'tolfun':1e-12}, 
        or the empty dict {}""" 
        res = {}
        if self.counteval > 0: 
            if self.counteval >= self.max_eval:
                res['evals'] = self.max_eval
            if self.ftarget is not None and len(self.fitvals) > 0 \
                    and self.fitvals[0] <= self.ftarget:
                res['ftarget'] = self.ftarget
            if max(self.D) > 1e7 * min(self.D):
                res['condition'] = 1e7
            if len(self.fitvals) > 1 \
                    and self.fitvals[-1] - self.fitvals[0] < self.tolfunhist:
                res['tolfunhist'] = self.tolfunhist
            if self.sigma * max(self.D) < self.tolx:
                # remark: max(D) >= max(diag(C))**0.5
                res['tolx'] = self.tolx
        return res

def BoundariesFunc(bounds, scaling_of_variables=None):
    _bounds = []
    if scaling_of_variables:
        for (a, b), s in zip(bounds, scaling_of_variables):
            _bounds.append( BoxContraintsTransformation(a/s, b/s) )
    else:
        for a, b in zip(bounds):
            _bounds.append( BoxContraintsTransformation(a, b) )
    return _bounds


def fmin(f, x0, sigma=1, **kwargs):
    """non-linear non-convex minimization procedure. 
    The functional interface to CMA-ES. 

    Parameters
    ==========
        `f`
            a function that takes as input a list of floats (like
            [3.0, 2.2, 1.1]) and returns a single float (a scalar).
            The objective is to find ``x`` with ``objectivefct(x)``
            to be as small as possible.
        `x0`
            list of numbers (like `[3,2,1.2]`), initial solution vector
        `sigma`
            float, initial step-size, standard deviation in any coordinate
        `tolfun`
            float, target function value, None for never
        `tolx`
            float, x-tolerance value, 1e-11 for never
        `tolfunhist`
            float, function-tolerance value, 1e-12 for never
        `max_eval`
            int or str, maximal number of function evaluations, a string
            is evaluated with N being the search space dimension
        `scaling_of_variables`
            list of numbers (like `[3,2,1.2]`), scaling vector
        `bounds`
            list of numbers (like `[(0,10), (0, 10), (0, 10)]`), boundary vector
        `verb_disp`
            int, display on console every verb_disp iteration, 0 for never
        `popsize`
            int, number of X per iteration,` 4 + int(3 * log(N))` for never
        `parallel`
            Boolean, provide `pop_size` X to the objective function, `False` for never
	    
           
    Returns
    =======
    ``return fbest, xbest
    """

    scaling_of_variables = kwargs.pop('scaling_of_variables', None)
    bounds = kwargs.pop('bounds', None)
    parallel = kwargs.pop('parallel', True)
    verb_disp = kwargs.pop('verb_disp', 1)
    if 'seed' in kwargs:
        random.seed(kwargs.pop('seed'))
    
    if bounds:
        bounds = BoundariesFunc(bounds, scaling_of_variables)
    
    def wrapper(func, scaling_of_variables=None, bounds=None):
        if bounds:
	   if scaling_of_variables:
	       if parallel:
	           def decorator(Xs):
		       Y = list()
		       for X in Xs:
	                   Y.append( [ s*f(x) for x, f, s in zip(X, bounds, scaling_of_variables) ] )
		       return func(Y)
	       else:
	           def decorator(X):
	               Y = [ s*f(x) for x, f, s in zip(X, bounds, scaling_of_variables) ]
		       return func(Y)
	   else:
	       if parallel:
	           def decorator(Xs):
		       Y = list()
		       for X in Xs:
	                   Y.append( [ f(x) for x, f in zip(X, bounds) ] )
		       return func(Y)
	       else:
	           def decorator(X):
	               Y = [ f(x) for x, f in zip(X, bounds) ]
		       return func(Y)
        else:
	   if scaling_of_variables:
	       if parallel:
	           def decorator(Xs):
		       Y = list()
		       for X in Xs:
	                   Y.append( [ s*x for x, s in zip(X, scaling_of_variables) ] )
		       return func(Y)
	       else:
	           def decorator(X):
	               Y = [ s*x for x, s in zip(X, scaling_of_variables) ]
		       return func(Y)
	   else:
	       decorator = func
        return decorator
    
    if scaling_of_variables:
	x0 = [ x/s for x, s in zip(x0, scaling_of_variables) ]
    if bounds:
	x0 = [ _f.inverse(x) for x, _f in zip(x0, bounds) ]
    if scaling_of_variables or bounds:
        f = wrapper(f, scaling_of_variables, bounds)

    opts = dict( max_eval = kwargs.get('max_eval', '1e3*N**2'),
                 ftarget = kwargs.get('tolfun', None),
		 tolx = kwargs.get('tolx', 1e-11),
		 tolfun = kwargs.get('tolfun', 1e-12),
		 tolfunhist = kwargs.get('tolfunhist', 1e-12),
		 popsize = kwargs.get('popsize', '4 + int(3 * log(N))'),
    )
    es = CMAES(x0, sigma, **opts)
    while not es.stop(): 
        X = es.ask()
	if parallel:
	    Y = f(X)
	else:
	    Y = [ f(x) for x in X]
        es.tell(X, Y)
        es.disp(verb_disp)
	
    if bounds:
        for i, f in enumerate(bounds):
            es.best.__dict__['x'][i] = f(es.best.__dict__['x'][i])
    if scaling_of_variables:
        for i, s in enumerate(scaling_of_variables):
            es.best.__dict__['x'][i] *= s
        
    print 'Termination by', 
    for k, v in  es.stop().iteritems():
        print k, '=', v,
    print
    #cma.pprint(es.best.__dict__)
    return es.best.__dict__['f'], list(es.best.__dict__['x'])
    
