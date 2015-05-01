from exceptions import StopIteration

class NewtonRaphson(object):
    """
        f : function for zeroing
	x0 : first x value for starting
	x1 : second x in case of secant method, default None
	epsfcn : used for calculate x1 in case of secant method, if x1 is not given, default 1e-4
	xtol : x tolerance for exit
	ftol : funtion tolerance for exit
	maxfev : maximum function evaluations for exit
	verbose : verbosity of iterations
	fprime : derivative function, default None
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
            raise StopIteration()
	    
	# update self.x
        self.x = self.x + self.damping*self.dx
        if self.bounds:
            if self.x<min(self.bounds):
                self.x = min(self.bounds)
            if self.x>max(self.bounds):
                self.x = max(self.bounds)
            
	
	# check stopping criteria
        if self.stopping_criteria():
            raise StopIteration()
	    
        return self.x, self.f(self.x)

    def run(self):
        if self.verbose:
            print 'x\t\tf(x)'
        i = 0
        for x, fx in self:
            if i>200:
                print 'not good'
                break
            if self.verbose:
                print '{x}\t{fx}'.format(x="%0.5e"%x, fx="%0.5e"%fx)
            i+=1
        if self.verbose:
            print '{x}\t{fx}'.format(x="%0.7e"%self.x, fx="%0.7e"%self.f(self.x))
        return self.x

    def stopping_criteria(self):
        if abs(self.dx) < self.xtol:
            return True
        if abs(self.f(self.x)) < self.ftol:
            return True
        if self.ifev > self.maxfev:
            return True

if __name__ == '__main__':

    from math import *
    def f(x):
        return cos(x)-x**3
 
    def f_prime(x):
        return -sin(x) - 3*x**2
                
    opt = NewtonRaphson(f=f, x0=0.5, verbose=True, bounds=(0.87, 0.88))
    x = opt.run()
    print 'result:', x


    
