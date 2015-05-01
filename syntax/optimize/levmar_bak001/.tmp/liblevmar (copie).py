# -*- coding: cp1252 -*-

from libarray  import *


class levmar(object):
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
    tau = 10.0e-3
    mu = 2.0
    eps_gradient = 1e-4
    eps_x_step = 1e-8
    iter = 0
    gradient = False
    x_step = False
    internal_mem = [None, None, None, None, None]
    external_criteria = lambda obj: False
    max_iter = 100
    function_call = 0
    stop = False
    dX = array([])
    
    def __init__(self):
        pass
        
    def config(self, **args):
        if 'damping' in args:  self.tau, self.mu = args['damping']
        if 'eps_gradient' in args:     self.eps_gradient = args['eps_gradient']
        if 'eps_x_step' in args:     self.eps_x_step = args['eps_x_step']
        if 'solve' in args:     self.solve = args['solve']
        if 'modarray' in args:
            modarray = __import__(args['modarray'])
            for comp in dir(modarray):
                globals()[comp]=getattr(modarray, comp)
                locals()[comp]=getattr(modarray, comp)
        if 'external_criteria' in  args: self.external_criteria = args['external_criteria']
        if 'max_iter' in  args: self.max_iter = args['max_iter']

    def update(self):
        # do update after step if step return True
        self.X += self.dX
        # apply new state
        self.X, self.A, self.g = self.get_memory(1).X, \
                                 self.get_memory(1).A, \
                                 self.get_memory(1).g
        return True

    def put_memory(self, slot):
        if slot>len(self.internal_mem) : raise 'memory slot is not available'
        self.internal_mem[slot] = __import__('copy').copy(self)
        return True

    def get_memory(self, slot):
        if slot>len(self.internal_mem) : raise 'memory slot is not available'
        return self.internal_mem[slot]

    def best_step(self):
        if self.get_memory(2)==None:
            self.put_memory(2)
        elif (norm(self.residual)<norm(self.get_memory(2).residual)):
            self.put_memory(2)
        else: pass
        return self.get_memory(2)

    def halting_criteria(self):
        ### halting criteria
        self.norm_residual = norm(self.residual)
        self.norm_dX = norm(self.dX)
        
        if 'residual' in dir(self):
            self.gradient = self.norm_residual<self.eps_gradient
        if 'dX' in dir(self):
            self.x_step = self.norm_dX<self.eps_x_step*(norm(self.X)+self.eps_x_step)

        return ( self.gradient \
                 or self.x_step \
                 or not(self.iter<self.max_iter) \
                 or self.external_criteria() )
           
    def init(self, X, residual, Jacobian):
        self.X = array([X.tolist()]).transpose()
        self.dX = self.X*0.0

        self.residual = residual(self)
        self.Jacobian = Jacobian(self)

        J_T = self.Jacobian.transpose()
        self.A = J_T.dot(self.Jacobian)
        self.g = -J_T.dot(self.residual)

        self.Fb = 0.5*self.residual.dot(self.residual)[0][0]

	self.halting_criteria()
        self.best_step()
        return True

    def step(self, residual, Jacobian):
        self.iter += 1

        # save init state
        self.put_memory(0)
        # solve delta X
        self.dX = self.solve(self.A + self.tau*self.A.diag(),self.g)

        # apply temporary solved iteration
        # to process errors and gain
        self.X += self.dX

        # call function and jacobian
        self.residual = residual(self)
        self.stop = self.halting_criteria()
        if self.stop : return False
        self.Jacobian = Jacobian(self)
        
        J_T = self.Jacobian.transpose()
        self.A = J_T.dot(self.Jacobian)
        self.g = -J_T.dot(self.residual)

        self.stop = self.halting_criteria()
        
        # save if best step
        self.best_step()

        self.Fn = 0.5*self.residual.dot(self.residual)[0][0]
        self.dL = (0.5*(self.dX.transpose().dot(self.tau*self.A.diag().dot(self.dX) + self.g)))[0][0]
        self.gain = (self.Fb-self.Fn)/self.dL

        # save new state
        self.put_memory(1)

        if self.stop : return False
        
        # restore init state
        self.X, self.A, self.g = self.get_memory(0).X, \
                                 self.get_memory(0).A, \
                                 self.get_memory(0).g

        # return self.gain>0 (i.e., if iteration is convergent or divergent)
        return (self.Fb>self.Fn and self.dL>0.0)

    def damping(self):
        """ calcul les parametres de damping pour chaque iteration
            convergente ou divergente
        """
        if self.stop : return False
        
        if self.Fb>self.Fn and self.dL>0.0:
            self.tau = self.tau*max(0.333333333, (1.0-(2.0*self.gain-1.0)**3)/2.0)
            self.mu = 2.0
            self.Fb = self.Fn
            return True
        else:
            self.tau = self.tau*self.mu
            self.mu = 2.0*self.mu
            return True
        return False


class LMA(levmar):
    def __init__(self, X, residual, Jacobian):
        self.X, self.func, self.Jac = X, residual, Jacobian
        wrap_gauss_pivot_total = __import__('libsolve').wrap_gauss_pivot_total
        self.config(eps_gradient=1e-6, eps_x_step=1e-6)
        self.config(damping=(1e-3, 2.0))
        self.config(solve=wrap_gauss_pivot_total)
        self.config(modarray='libarray')
	self.echo = lambda obj, step: None 
	self.out = lambda obj: (obj.X, obj.norm_residual, obj.norm_dX)   
           
    def run(self):
        self.echo(self, 0)
        self.init(self.X, self.func, self.Jac)
        self.echo(self, 1)
        while not(self.stop):
            if self.step(self.func, self.Jac):
                self.update()
            self.damping()
            self.echo(self, 2)
        self.bstep= self.best_step()
        self.echo(self, 3)
        return self.out(self.bstep)


       
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
        print '%7d'%obj.bstep.iter,
        for v in obj.bstep.X.transpose()[0].tolist():
            print '%11.3e'% v, 
        print '%11.3e'%obj.bstep.norm_residual,
        print
        print '-- function call : %d -----------' % (obj.function_call)
        return 
            


