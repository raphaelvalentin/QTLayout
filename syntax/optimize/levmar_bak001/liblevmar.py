# -*- coding: cp1252 -*-
from math import *
from exceptions import *
from libarray  import *

         
def norm(vec):
    return sqrt(sum([x*x for x in vec]))


class levmar(list):
    """
    * Classe generaliste formant les constituants
    de l_algorithme d_optimisation de Levenberg - Marquart.
    * Necessite une classe 'array' (numpy/libarray) et un solveur
      lineaire nomme 'solve'.
    * Parametres config :
     - damping : parametres de damping par defaut (1e-3, 2.0).
                 Le couple (100e-3, 2.0) peut augmenter la convergence
                 en contrepartie de la vitesse de convergence.
     - eps_gradient : erreur absolue minimal de la fonction f-f0 pour le critere d_arret
     - eps_x_step : distance minimale entre deux steps pour le critere d_arret
     - dp     = fractional increment of 'p' for numerical derivatives
    * additional config :
     - external_criteria : critere externe d_arret, parametre : fonction(obj)/lambda obj,
                          retourne un booleen.
     * input:
        X : array([...]) vecteur initial X
        residual : fonction d_erreur (residual) f0-f(...), argument : objet levmar,
            retourne : array (1D)
        Jacobian : fonction jacobienne, objet levmar attendu en argument
            retourne : array (2D)
    * exemple d_implentation de l_algorithme de LM:
    
        def residual(objet):
            X = objet.X
            ...
            return array([,...])

        def Jacobian(objet):
            X = objet.X
            f0 = objet.f
            ...
            return array([[,...],...])
        
        optim = levmar()
        
        optim.config(eps_gradient=1e-4, eps_x_step=1e-4, max_iter=100)
        optim.config(damping=(1.0e-3, 2.0))
        optim.config(solve=gauss_pivot_total)
        optim.config(modarray='numpy')
        
        optim.init(X, residual, Jacobian)
        while halting_criteria() is False:
            if optim.step(f, J):
                optim.update()
            optim.damping()
            print optim.iter, optim.X

    """
    
    ifev = 0
    fev = []
    def wrapf(self, func):
    # count + memoize the evals
        def call(*x):
            for _x, fx in self.fev:
                if _x == x:
                    return fx
            fx = func(*x)
            self.ifev += 1
            fev.append((x, fx))
            return fx
        return call

    gradient = 1e300
    xstep = 1e300
    stop = False
    
    def __init__(self, func, x0, Dfun=None, ftol=1.0e-04, xtol=1.0e-08, maxfev=100, epsfcn=1e-6, damping=(10e-3, 2.0)):
        self.func = self.wrapf(self, func)
        self.x0 = x0
        if Dfun:
            self.Dfun = Dfun
        self.ftol = ftol
        self.xtol = xtol
        self.maxfev = maxfev
        self.epsfcn = epsfcn
        self.damping = damping
        self.k = 0

        self.tau, self.nu = self.damping
        self.X = array([float(xi) for xi in x0])
        
    def update(self):
        # do update after step if step return True
        self.X += self.dX

    def Dfun(self, x0):
        """ numerical jacobian matrix
        """
        J = array([])
        f0 = self.func(x0)
        if '__iter__' in dir(self.epsfcn):
            for i in xrange(len(x)):
                X_dx = [xi for xi in x]
                X_dx[i] += self.epsfcn[i]
                J.append(-(self.func(X_dx) - f0)/self.epsfcn[i])
        else:
            for i in xrange(len(x)):
                X_dx = [xi for xi in x]
                X_dx[i] += self.epsfcn
                J.append(-(self.func(X_dx) - f0)/self.epsfcn)
        return J.T

    def __enter__(self):
        # init
        f = self.func(self.X)
        J = self.Dfun(self.X)
        self.A = J.T.dot(J)
        self.g = -J.T.dot(f)
        self.mu = self.tau*self.A.diag()
        return self

    def best_step(self):
        if self.get_memory(2)==None:
            self.put_memory(2)
        elif (norm(self.residual)<norm(self.get_memory(2).residual)):
            self.put_memory(2)
        else: pass
        return self.get_memory(2)

    def stopping_criteria(self):
        ### stopping criteria
        self.gradient = norm(self.g)
        if self.iter>0:
            self.xstep = norm(self.dX)

        if (self.gradient<self.ftol) \
           or (self.xstep<self.xtol*(norm(self.X)+self.xtol)) \
           or (self.ifev>self.maxfev):
           raise StopIteration()
           

    def step(self):
        # save init state
        self.put_memory(0)
        # solve delta X
        I = identity(len(A))
        self.dX = self.solve(self.A + dot(self.mu, I), self.g)

        # apply temporary solved iteration
        # to process errors and gain
        self.X += self.dX

        # call function and jacobian
        self.residual = residual(self)
        self.stop = self.halting_criteria()
        if self.stop : return False
        self.Jacobian = Jacobian(self)
        
        self.A = self.Jacobian.T.dot(self.Jacobian)
        self.g = -self.Jacobian.T.dot(self.residual)
        self.stop = self.halting_criteria()
        
        # save if best step
        self.best_step()

        self.Fn = 0.5*self.residual.dot(self.residual)[0][0]
        self.dL = (0.5*(self.dX.transpose().dot(self.tau*self.A.diag().dot(self.dX) + self.g)))[0][0]
        self.gain = (self.F0-self.Fn)/self.dL

        # save new state
        self.put_memory(1)

        if self.stop : return False
        
        # restore init state
        self.X, self.A, self.g = self.get_memory(0).X, \
                                 self.get_memory(0).A, \
                                 self.get_memory(0).g

        # return self.gain>0 (i.e., if iteration is convergent or divergent)
        return (self.F0>self.Fn and self.dL>0.0)
        
    def F(self, X):
        fx = self.func(X)
        return 0.5*dot(fx.T, fx)[0][0]    
        

    def damping(self):
        """ calcul les parametres de damping pour chaque iteration
            convergente ou divergente
        """
        if self.stop : return False
        
        if self.F0>self.Fn and self.dL>0.0:
            self.mu = self.mu*max(0.333333333, (1.0-(2.0*self.gain-1.0)**3)/2.0)
            self.mu = 2.0
            self.F0 = self.Fn
            return True
        else:
            self.tau = self.tau*self.mu
            self.mu = 2.0*self.mu
            return True
        return False

def LMA(*args, **kwargs):
    with levmar(*args, **kwargs) as opt:
        for step in opt:
            if step:
                opt.update()
            opt.damping()
    return opt.X, opt.norm_residual
                    
def echo(obj, step):
    if step==0:
        print '%7.8s'%'iter',
        for i in xrange(len(obj.X)):
            print '%12s'%'X%d'% i, 
        print '%11s'%'|residual|',
        print
        print '----------------------------------'
        return 
    if step==1 or step==2 :
        print '%7d'%obj.iter,
        for v in obj.X.transpose()[0].tolist():
            print '%11.3e'% v, 
        print '%11.3e'%obj.norm_residual,
        print
        return 
    if step==3:
        print '-- final ----------------------------'
        print '%12s'%obj.bstep.iter,
        for v in obj.bstep.X.transpose()[0].tolist():
            print '%11.3e'% v, 
        print '%11.3e'%obj.bstep.norm_residual,
        print
        print '-- function call : %d -----------' % (obj.function_call)
        return 
            


