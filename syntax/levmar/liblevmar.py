# -*- coding: cp1252 -*-
from math import *
from exceptions import *
#from numpy  import *
#from numpy import zeros, identity
#from numpy.linalg import solve, LinAlgError

from libarray import *
from libarray import zeros, identity
from linalg import solve, norm, LinAlgError, det


class LevMarError(Exception):
    pass

class levmar(list):
    """
    * Classe generaliste formant les constituants
    de l_algorithme d_optimisation de Levenberg - Marquart.
    * Necessite une classe 'array' (numpy/libarray) et un solveur
      lineaire nomme 'solve'.
    * exemple d_implentation de l_algorithme de LM:
    
        def f(x):
            [x0, x1...] = x
            ...
            return array([,...])
        x0 = [0., 0.]
        with levmar(f, x0) as opt:
            for x, residual in opt:
                pass
    """
    def wfunc(self, func):
    # count + memoize the evals
        self.ifev = 0
        self.fev = []
        def call(x, *args):
            for _x, _args, fx in self.fev:
                if reduce(lambda x, y: x and y, [v1==v2 for v1, v2 in zip(_x, x)]) and _args==args:
                    return fx
            if self.ifev<self.maxfev:
                fx = array(func(x, *args))
                self.ifev += 1
                self.fev.append((x, args, fx))
                return fx
            else:
                raise StopIteration()
        return call

    
    def __init__(self, func, x0, Dfun=None, bounds=[], gtol=1.0e-06, xtol=1.0e-06, maxfev=1000, epsfcn=1e-6, damping=(10e-3, 2.0), maxiter=100):
        self.func = self.wfunc(func)
        self.x0 = x0
        if Dfun:
            self.Dfun = Dfun
        self.gtol = gtol
        self.xtol = xtol
        self.maxfev = maxfev
        self.epsfcn = epsfcn
        self.damping = damping
        self.iter = 0
	self.maxiter = maxiter

        self.tau, self.nu = self.damping
        self.X = array([float(x) for x in x0])
	
	self.suppressExceptions = False
	self.traceback = ''
	if len(bounds)<>0 and len(bounds)<>len(self.x0):
	    raise LevMarError('The dimensions of bounds are not correct.')
	self.bounds = bounds
	

    def Dfun(self, x0):
        """ numerical jacobian matrix
        """
        f0 = self.func(x0)
        shape = (len(f0), len(x0))
        J = zeros(shape)
        if isinstance(self.epsfcn, (list, tuple)):
            for i in xrange(shape[1]):
                x1 = [x for x in x0]
                x1[i] += self.epsfcn[i]
                f1 = self.func(x1)
                for j in xrange(shape[0]):
                    J[j][i] = (f1[j]-f0[j])/self.epsfcn[i]
        else:
            for i in xrange(shape[1]):
                x1 = [x for x in x0]
                x1[i] += self.epsfcn
                f1 = self.func(x1)
                for j in xrange(shape[0]):
                    J[j][i] = (f1[j]-f0[j])/self.epsfcn
        return J

    def __enter__(self):
        self.f = self.func(self.X)
        f = array([self.f]).T
        self.J = self.Dfun(self.X)
        self.A = self.J.T.dot(self.J)
        self.g = -self.J.T.dot(f)
        I = identity(len(self.A))
        self.mu = self.tau*diag(self.A)*I        
        self.F0 = 0.5*f.T.dot(f)[0][0]
        self.append((self.iter, self.X, norm(self.f)))
	return self

    def __exit__(self, type, value, traceback):
        if isinstance(value, LevMarError): 
	    self.traceback = (type, value, traceback)
        if isinstance(value, LinAlgError): 
	    self.traceback = (type, value, traceback)
        return isinstance(value, LevMarError) or isinstance(value, LinAlgError) 

    def checkJacobian(self):
        shape = (len(self.J), len(self.J[0]))
        isJnull = [0.0]*shape[1]
        for j in xrange(shape[1]):
            for i in xrange(shape[0]):
                isJnull[j] += self.J[i][j]**2
        return tuple(sqrt(isJnull))

    def __iter__(self):
        return self

    def next(self):

	isJnull = self.checkJacobian()
	for i, f in enumerate(self.f):
	    if f == float('nan'):
	        raise LevMarError('The point %d is a Nan. The optimizer has been stopped.'%i)
	    
    
        # stopping criteria
        if norm(self.g.T[0])<self.gtol:
            raise StopIteration('Magnitude of gradient smaller than the \'gtol\' tolerance.')
        # increment iter
        self.iter += 1
        if self.iter>self.maxiter:
            raise StopIteration('Number of iterations exceeded \'maxiter\' or number of function evaluations exceeded \'maxfev\'.')
     
        
        # save previous state
        self.X0, self.A0, self.g0, self.f0 = self.X.copy(), self.A.copy(), self.g.copy(), self.f.copy(), 
	
        # compute dX
        I = identity(len(self.A))
        self.dX = array(solve(self.A + self.mu, self.g.T[0], verbose=False))
	
        # stopping criteria
        #if norm(self.dX)<self.xtol*(norm(self.X)+self.xtol):
        #    raise StopIteration()
	if not(False in [a<b for a, b in zip(abs(self.dX), self.xtol*(abs(self.X)+self.xtol))]):
            raise StopIteration('Change in x smaller than the \'xtol\' tolerance.')

        self.X = self.X + self.dX
        
        # compute f using the new X
        self.f = self.func(self.X)
        f = array([self.f]).T
        
        self.Fn = 0.5*f.T.dot(f)[0][0]
        
        dX = array([self.dX]).T
        self.dL = (0.5*(dX.T.dot(self.mu.dot(dX) + self.g)))[0][0]
	self.dF = self.F0-self.Fn
	
	if len(self.bounds):
	    isStepAcceptable = not(False in [min(bound)<=x<=max(bound) for x, bound in zip(self.X, self.bounds)])
	else:
	    isStepAcceptable = True
        
        # if step acceptable
        if ((self.dF>0 and self.dL>0) or (self.dF<0 and self.dL<0)) and isStepAcceptable:
            # compute jacobian, A and g
            self.J = self.Dfun(self.X)
            self.A = self.J.T.dot(self.J)
            self.g = -self.J.T.dot(f)

            # damp mu and update F0 parameter
            self.mu = self.mu*max(0.333333333, (1.0-(2.0*(self.dF/self.dL)-1.0)**3)/2.0)
            self.nu = self.damping[1]
            self.F0 = self.Fn
         
        else:
            # restore
            self.X, self.A, self.g, self.f = self.X0.copy(), self.A0.copy(), self.g0.copy(), self.f0.copy(), 
            # damp mu
            self.mu = self.mu*self.nu
            self.nu = 2.0*self.nu
                
        # save and return X and residual
        residual = norm(self.f)
        self.append((self.iter, self.X, residual))
        return self.iter, self.X, residual
       

def fmin(func, x0, Dfun=None, gtol=1.0e-06, xtol=1.0e-04, maxfev=1000, maxiter=50, epsfcn=1e-6, damping=(10e-3, 2.0), verbose=True):
        with levmar(func, x0, Dfun=Dfun, gtol=gtol, xtol=xtol, maxfev=maxfev, maxiter=maxiter, epsfcn=epsfcn, damping=damping) as opt:
            if verbose:
                print 'Start Levenberg Marquart Optimizer...'
	        print '{step:>7}{x}{residual:>13}'.format(step='step', x='{:>13}'*len(x0), residual='residual').format(*('X[%d]'%i for i in xrange(len(x0))))
	        print '{step:>7}{x}{residual:>13.3e}'.format(step=0, x='{:>13.3e}'*len(x0), residual=norm(opt.f)).format(*x0)
	    while 1:
	        iter=0
	        try:
		    iter, X, residual = opt.next()
                    if verbose:
	                print '{step:>7}{x}{residual:>13.3e}'.format(step=iter, x='{:>13.3e}'*len(X), residual=residual).format(*X)
		except StopIteration, exit_message:
		    if verbose:
		        print '---------'
		        print 'Optimization terminated successfully.'
		        print '         Exit Message: %s'%(exit_message)
		        print '         Iterations: %d'%(iter)
			print '         Function evaluations: %d'%(opt.ifev)
		    break
        if verbose:
            print '         Best step:'
	    if len(opt):
	        iter, X, residual  = min(list(opt[i] for i in xrange(len(opt))), key=lambda step: step[2])
                print '{step:>7}{x}{residual:>13.3e}'.format(step=iter, x='{:>13.3e}'*len(X), residual=residual).format(*X)
	    else:
	        iter = 0; residual = 0.0; X=x0
                print '{step:>7}{x}{residual:>13}'.format(step=0, x='{:>13.3e}'*len(x0), residual=norm(opt.f)).format(*x0)
	if opt.traceback:
	    print 'LevMarError:', opt.traceback[1]
	if verbose:
	    print 
        return iter, X, residual


def fit(func, x0, xmeas, ymeas, gtol=1.0e-08, xtol=1.0e-06, maxfev=100, epsfcn=1e-8, damping=(10e-3, 2.0), verbose=True):
    def f(X, xmeas, ymeas):
        return array([(ymeas-func(x, *X))/ymeas for x in xmeas])
    return fmin(f, x0, Dfun=None, gtol=gtol, xtol=xtol, maxfev=maxfev, epsfcn=epsfcn, damping=damping, verbose=verbose)
    
    
        
    
    


if __name__ == '__main__':

    from random import random
    noise = [(random()*2-1)*1e-9 for i in xrange(31)]
    xmeas = array([-0.5 ,-0.45 ,-0.4 ,-0.35 ,-0.3 ,-0.25 ,-0.2 ,-0.15 ,-0.1 ,-0.05 ,0.0 ,0.05 ,0.1 ,
                    0.15 ,0.2 ,0.25 ,0.3 ,0.35 ,0.4 ,0.45 ,0.5 ,0.55 ,0.6 ,0.65 ,0.7 ,0.75 ,0.8 ,0.85 ,0.9 ,0.95 ,1.0])
    def f(x, u0, ua, ub, vth):
        tox = 1e-1
        return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
    ymeas = array([f(x, 300, 0.03, 0.04, 0.2) for x in xmeas])
    ymeas += noise

    def f(X):
        [u0, ua, ub, vth] = X
        tox = 1e-1 
        y = u0/(1.0+(ua/tox)*(xmeas-vth)+(ub/tox/tox)*(xmeas-vth)**2)
        return ((ymeas-y)/ymeas)**2 + ((u0-100)/100.)**2
        
    X =array([500., 0.1, 0.00, 0.25])
    nfev, X, residual = fmin(f, X, epsfcn=1e-4, gtol=1.0e-06, xtol=1.0e-06, )
    
    exit()
    def f(x, u0, ua, ub, vth):
        tox = 1e-1
        return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
    
    #xmeas = array([-0.5 ,-0.45 ,-0.4 ,-0.35 ,-0.3 ,-0.25 ,-0.2 ,-0.15 ,-0.1 ,-0.05 ,0.0 ,0.05 ,0.1 ,0.15 ,0.2 ,0.25 ,0.3 ,0.35 ,0.4 ,0.45 ,0.5 ,0.55 ,0.6 ,0.65 ,0.7 ,0.75 ,0.8 ,0.85 ,0.9 ,0.95 ,1.0])
    #ymeas = array([108.7 ,120.0 ,132.74 ,147.06 ,163.04 ,180.72 ,200.0 ,220.59 ,241.94 ,263.16 ,283.02 ,300.0 ,312.5 ,319.15 ,319.15 ,312.5 ,300.0 ,283.02 ,263.16 ,241.94 ,220.59 ,200.0 ,180.72 ,163.04 ,147.06 ,132.74 ,120.0 ,108.7 ,98.68 ,89.82 ,81.97])
    #nfev, X, residual = fit(f, X, xmeas, ymeas)

