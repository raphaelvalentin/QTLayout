from exceptions import StopIteration

class NewtonRaphson(object):
    """
Newton Raphson 1D Local Optimizer
\tf : function for zeroing
\tx0 : first x value for starting
\tx1 : second x in case of secant method, default None
\tepsfcn : used for calculate x1 in case of secant method, if x1 is not given, default 1e-4
\txtol : x tolerance for exit
\tftol : funtion tolerance for exit
\tmaxfev : maximum function evaluations for exit
\tverbose : verbosity of iterations
\tfprime : derivative function, default None
    """

    # wrap function to add memoize and counter properties
    ifev = 0
    cache = []
    def wrap(self, f):
        def call(*x):
            for xi, fxi in self.cache:
                if x == xi:
                    return fxi
            fx = f(*x)
            self.ifev += 1
            self.cache.append((x, fx))
            return fx
        return call

    def __init__(self, f, x0, **kwargs):
        self.kwargs = kwargs
        self.f = self.wrap(f)
        self.x = x0
        self.cache = []
        self.epsfcn = kwargs.get('epsfcn', 1e-4)
        self.xtol = kwargs.get('xtol', 1e-6)
        self.ftol = kwargs.get('ftol', 1e-6)
        self.maxfev = kwargs.get('maxfev', 100)
        self.verbose = kwargs.get('verbose', True)
        self.bounds = kwargs.get('bounds', tuple())
        self.damping = abs(kwargs.get('damping', 1.0))
        self.iter = 0
        if 'fprime' in kwargs:
            self.fprime = kwargs['fprime']
        else:
            if 'x1' in kwargs:
                self.x0 = self.x
                self.x =  kwargs['x1']
            else:
                if self.x >= 0:
                    self.x0 = self.x*(1.0 - self.epsfcn) - self.epsfcn
                else:
                    self.x0 = self.x*(1.0 - self.epsfcn) + self.epsfcn

    def __iter__(self):
        return self
    
    def next(self):
        # evaluate the progression
        try:
            if 'fprim' in dir(self):
                # Newton method Step
                self.dx = -self.f(self.x)/self.fprime(self.x)
            else:
                # Secant method Step
                self.dx = -self.f(self.x)*(self.x - self.x0)/(self.f(self.x) - self.f(self.x0))
                self.x0 = self.x
        except ZeroDivisionError:
            raise StopIteration('ZeroDivisionError')
	    
	# update self.x
        self.x = self.x + self.damping*self.dx
        if self.bounds:
            if self.x<min(self.bounds):
                self.x = min(self.bounds)
            if self.x>max(self.bounds):
                self.x = max(self.bounds)
            
	
	# check stopping criteria
        if self.stopping_criteria():
            raise StopIteration('')

        self.iter += 1
        return self.iter, self.x, self.f(self.x)

    def run(self):
        if self.verbose:
            print 'Start Newton-Raphson Optimizer...'
            print '{step:>7}{x:>13}{fx:>13}'.format(step='step', x='x', fx='f(x)')
            print '{step:>7}{x:>13.3e}{fx:>13.3e}'.format(step=0, x=self.x, fx=self.f(self.x))

        while 1:
            iter=0
            try:
                iter, x, fx = self.next()
                if self.verbose:
                    print '{step:>7}{x:>13.3e}{fx:>13.3e}'.format(step=iter, x=x, fx=fx)
            except StopIteration, exit_message:
                if self.verbose:
                    print '---------'
                    if exit_message <> 'ZeroDivisionError':
                        print 'Optimization terminated successfully.'
                    else:
                        print 'Optimization terminated accidently!'
                    print '         Exit Message: %s'%(exit_message)
                    print '         Iterations: %d'%(self.iter)
                    print '         Function evaluations: %d'%(self.ifev)
                break
        if self.verbose:
            print '         Best step:'
            print '{step:>7}{x:>13.3e}{fx:>13.3e}'.format(step=self.iter, x=self.x, fx=self.f(self.x))
        return self.x

    def stopping_criteria(self):
        if abs(self.dx) < self.xtol:
            raise StopIteration('Change in x smaller than \'xtol\' tolerance.')
        if abs(self.f(self.x)) < self.ftol:
            raise StopIteration('Magnitude of \'f(x)\' smaller than \'ftol\' tolerance.')
        if self.ifev > self.maxfev:
            raise StopIteration('Number of function evaluations exceeded \'maxfev\'.')


def fmin(*args, **kwargs):
    """
Newton Raphson 1D Local Optimizer
\tf : function for zeroing
\tx0 : first x value for starting
\tx1 : second x in case of secant method, default None
\tepsfcn : used for calculate x1 in case of secant method, if x1 is not given, default 1e-4
\txtol : x tolerance for exit
\tftol : funtion tolerance for exit
\tmaxfev : maximum function evaluations for exit
\tverbose : verbosity of iterations
\tfprime : derivative function, default None
    """
    return NewtonRaphson(*args, **kwargs).run()


if __name__ == '__main__':

    from math import *
    def f(x):
        return cos(x)-x**3
 
    def f_prime(x):
        return -sin(x) - 3*x**2
                
    #opt = NewtonRaphson(f=f, x0=0.5, verbose=True, bounds=(0.0, 1.0))
    opt = fmin(f=f, x0=0.5, verbose=True, bounds=(0.0, 1.0))
    #x = opt.run()
    #print 'result:', x


    
