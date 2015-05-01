# -*- coding: cp1252 -*-
from math import *
from exceptions import *

from libarray import *
from libarray import zeros, identity, nan
from linalg import solve, norm, LinAlgError, det, chol
from optimize.transformations import BoxContraintsTransformation, ScalingTransformation


__all__ = ['levmar', 'fmin', 'LevMarWarning']

class LevMarError(Exception):
    pass

class LevMarWarning:
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return str(self.message)    

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
    def wrapper_func(self, func):
    # count + memoize the evals
        self.ifev = 0
        self.fev = []
        def call(x, *args):
            for _x, _args, fx in self.fev:
                if reduce(lambda x, y: x and y, [v1==v2 for v1, v2 in zip(_x, x)]) and _args==args:
                    return fx
            if self.ifev<self.maxfev:
	        if len(self.bounds)<>0:
		    x = type(x)([f(x) for x, f in zip(x, self.box_contraints_transformation)])
	        if len(self.scaling_of_variables)<>0:
		    x = type(x)([f.inverse(x) for x, f in zip(x, self.scaling_transformation)])
                fx = func(x, *args)
                self.ifev += 1
                self.fev.append((x, args, fx))
                return fx
            else:
                raise StopIteration()
        return call

    
    def __init__(self, func, x0, Dfun=None, bounds=[], scaling_of_variables=[], gtol=1.0e-06, xtol=1.0e-06, maxfev=1000, epsfcn=1e-6, damping=(1e-3, 2.0), maxiter=100):

        self.func = func
        self.x0 = x0
        if Dfun:
            self.Dfun = Dfun
	self.bounds = bounds if bounds<>None else []
	self.scaling_of_variables = scaling_of_variables if scaling_of_variables<>None else []
        self.gtol = gtol
        self.xtol = xtol
        self.maxfev = maxfev
        self.epsfcn = epsfcn
        self.damping = damping
	self.maxiter = maxiter
        self.tau, self.nu = self.damping
	
        self.iter = 0
	self.suppressExceptions = False
	self.traceback = ''
	
	if len(self.bounds)<>0 and len(self.bounds)<>len(self.x0):
	    raise LevMarError('The dimensions of bounds are not correct.')
	if len(self.scaling_of_variables)<>0:
	    if len(self.scaling_of_variables)<>len(x0):
	        raise LevMarError('The dimensions of scaling_of_variables are not correct.')
	    
	self.scaling_transformation = [ScalingTransformation(s) for s in self.scaling_of_variables]
	self.box_contraints_transformation = []
	
	if len(self.scaling_of_variables)<>0:
	    for b, fsov in zip(self.bounds,self.scaling_transformation):
	        self.box_contraints_transformation.append( BoxContraintsTransformation( fsov(bi) for bi in b ) )
	else:
	    for b in self.bounds:
	        self.box_contraints_transformation.append( BoxContraintsTransformation(b) )
	
	if len(self.scaling_of_variables)<>0:
	    if len(self.bounds)<>0:
	        self.X = array([b.inverse(f(float(x))) for x, f, b in zip(self.x0, self.scaling_transformation, self.box_contraints_transformation)])
	    else:
	        self.X = array([f(float(x)) for x, f in zip(self.x0, self.scaling_transformation)])
	else:
	    if len(self.bounds)<>0:
	        self.X = array([b.inverse(float(x)) for x, b in zip(self.x0, self.box_contraints_transformation)])
	    else:
	        self.X = array([float(x) for x in self.x0])

        self.func = self.wrapper_func(func)
	

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
        self.mu = self.tau
        self.F0 = 0.5*f.T.dot(f)[0][0]
	
	X = self.X
	if len(self.bounds)<>0:
	    X = type(X)([f(x) for x, f in zip(X, self.box_contraints_transformation)])
	if len(self.scaling_of_variables)<>0:
	    X = [f.inverse(x) for x, f in zip(X, self.scaling_transformation)]
	self.append((self.iter, X, norm(self.f)))
	return self

    def __exit__(self, type, value, traceback):
        if isinstance(value, LevMarError): 
	    self.traceback = (type, value, traceback)
        if isinstance(value, LinAlgError): 
	    self.traceback = (type, value, traceback)
        return isinstance(value, LevMarError) or isinstance(value, LinAlgError) 

    def checkJacobian(self):
        shapeJ = shape(self.J)
        isJnull = [0.0]*shapeJ[1]
        for j in xrange(shapeJ[1]):
            for i in xrange(shapeJ[0]):
                isJnull[j] += self.J[i][j]**2
        return tuple(sqrt(isJnull))

    def __iter__(self):
        return self

    def next(self):

	isJnull = self.checkJacobian()
	if 0 in isJnull:
	    print LevMarWarning('The Jacobian is rank deficient.')
	
	if nan in self.f:
            raise LevMarError('One point of \'f\' is a Nan. The optimizer has been stopped.')

        # stopping criteria
        if norm(self.g)<self.gtol:
            raise StopIteration('Magnitude of gradient smaller than the \'gtol\' tolerance.')
	
        # increment iter
        self.iter += 1
        if self.iter>self.maxiter:
            raise StopIteration('Number of iterations exceeded \'maxiter\' or number of function evaluations exceeded \'maxfev\'.')
     
        # save previous state
        self.X0 = self.X.copy()
	self.A0 = self.A.copy()
	self.g0 = self.g.copy()
	self.f0 = self.f.copy()
	
        # compute dX
        I = identity(len(self.A))
	A = self.A
	mu = self.mu
	g = self.g
	
	def geth(A, g, mu):
	    mA = max(abs(A.flatten()))
	    for i in xrange(5):
	        try:
		    I = identity(len(A))
		    L = chol(A + mu*diag(A)*I)
		    h = L.solve(g)
		    return array(h), mu
		except LinAlgError:
		    mu = max(10.0*mu, 1e-15*mA)
            raise LevMarError('The matrix A + mu*diag(A) is not positive. The optimizer has been stopped.')
	    
        #self.dX = array(solve(A + mu*diag(A)*I, g.T[0]))
	self.dX, self.mu = geth(A, g.T[0], mu)
	
        # stopping criteria
        if norm(self.dX)<self.xtol*(norm(self.X)+self.xtol):
            raise StopIteration('Change in x smaller than the \'xtol\' tolerance.')

        self.X = self.X + self.dX
        
        # compute f using the new X
        self.f = self.func(self.X)
        f = array([self.f]).T
        dX = array([self.dX]).T
	
        self.Fn = 0.5*f.T.dot(f)[0][0]
        
	# Using the state A, mu, g, h and considering a Taylor expansion model L, the trust region dL = L(h)-L(0) is calculated.
	# then dF/dL is the gain of the trust region. dL is positive
        self.dL = 0.5*dX.T.dot((self.mu*diag(self.A)*I).dot(dX) + self.g)[0][0]
	self.dF = self.F0-self.Fn
	
        # if step acceptable, dF is sometimes null due to machine accuracy
        if self.dF>0 and self.dL>0 :
            # compute jacobian, A and g
            self.J = self.Dfun(self.X)
            self.A = self.J.T.dot(self.J)
            self.g = -self.J.T.dot(f)

            # damp mu and update F0 parameter
            self.mu = self.mu*max(1.0/3.0, (1.0-(2.0*(self.dF/self.dL)-1.0)**3)/2.0)
            self.nu = self.damping[1]
            self.F0 = self.Fn
	    
	    if norm(self.g0-self.g) < self.gtol:
                raise StopIteration('Change in g smaller than the \'gtol\' tolerance.')

        else:
            # restore
            self.X = self.X0.copy()
	    self.A = self.A0.copy()
	    self.g = self.g0.copy()
	    self.f = self.f0.copy()
            # damp mu
            self.mu = self.mu*self.nu
            self.nu = 2.0*self.nu


        # save and return X and residual
        residual = norm(self.f)
	
	X = self.X
	if len(self.bounds)<>0:
	    X = type(X)([f(x) for x, f in zip(X, self.box_contraints_transformation)])
	if len(self.scaling_of_variables)<>0:
	    X = [f.inverse(x) for x, f in zip(X, self.scaling_transformation)]

	self.append((self.iter, X, residual))
	return self.iter, X, residual
       




def fmin(func, x0, Dfun=None, bounds=[], scaling_of_variables=[], gtol=1.0e-06, xtol=1.0e-04, maxfev=1000, maxiter=100, epsfcn=1e-6, damping=(1000e-3, 2.0), verbose=True):
        with levmar(func, x0, Dfun=Dfun, bounds=bounds, scaling_of_variables=scaling_of_variables, gtol=gtol, xtol=xtol, maxfev=maxfev, maxiter=maxiter, epsfcn=epsfcn, damping=damping) as opt:
            if verbose:
                print 'Start Levenberg Marquart Optimizer...'
	        print '{step:>7}{x}{residual:>13}'.format(step='step', x='{:>13}'*len(x0), residual='residual').format(*('X[%d]'%i for i in xrange(len(x0))))
	        print '{step:>7}{x}{residual:>13.3e}'.format(step=0, x='{:>13.3e}'*len(x0), residual=norm(opt.f)).format(*x0)
            iter = 0		
	    while 1:
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


  
    
    


if __name__ == '__main__':

    from random import gauss, seed
    from function import linspace
    #seed(1234)
    noise = array([gauss(0, 1) for i in xrange(101)])
    xmeas = linspace(0.2, 1.2, 101)
    def f(x, u0, ua, ub, vth):
        tox = 1e-1
        mu = u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
	return mu*(x-vth)
    ymeas = f(xmeas, 300, 0.1, 0.1, 0.2) 
    ymeas += noise

    def f(X):
        [u0, ua, ub, vth] = X
        tox = 1e-1 
        mu = u0/(1.0+(ua/tox)*(xmeas-vth)+(ub/tox/tox)*(xmeas-vth)**2)
	y = mu*(xmeas-vth)
        return (ymeas-y)
        
    X =array([500., 0.00, 0.00, 0.5])
    bounds = [(50, 600),(0, 1),(0 ,1),(-1, 1)]
    #bounds = []
    scaling_of_variables=array([300., 0.1, 0.1, 0.1])
    #scaling_of_variables=[]
    nfev, X, residual = fmin(f, X, epsfcn=1e-4, gtol=1.0e-6, xtol=1.0e-6, scaling_of_variables=scaling_of_variables, bounds=bounds, maxiter=200)
    
    def f(x, u0, ua, ub, vth):
        tox = 1e-1
        mu = u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
	return mu*(x-vth)
    
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(xmeas, ymeas, color='black', label='meas')
    ysim = f(xmeas, *X)
    plt.plot(xmeas,ysim, color='red', label='opt')
    plt.legend(loc=0)
    plt.show()

    
