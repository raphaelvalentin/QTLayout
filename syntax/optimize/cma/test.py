#https://www.lri.fr/~hansen/cma.py
import matplotlib.pyplot as plt
from numpy import *
from optimize import cma

#plt.plot([1,2,3], [2,8,1])
#plt.show()

from optimize import cma


from random import random
noise = array([(random()*2-1)*1e-20 for i in xrange(31)])
xmeas = array([-0.5 ,-0.45 ,-0.4 ,-0.35 ,-0.3 ,-0.25 ,-0.2 ,-0.15 ,-0.1 ,-0.05 ,0.0 ,0.05 ,0.1 ,
                0.15 ,0.2 ,0.25 ,0.3 ,0.35 ,0.4 ,0.45 ,0.5 ,0.55 ,0.6 ,0.65 ,0.7 ,0.75 ,0.8 ,0.85 ,0.9 ,0.95 ,1.0])
def f(x, u0, ua, ub, vth):
    tox = 1e-1
    return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
ymeas = array([f(x, 300, 0.1, 0.1, 0.2) for x in xmeas])
ymeas += noise

def f(X):
    [u0, ua, ub, vth] = X
    tox = 1e-1 
    y = u0/(1.0+(ua/tox)*(xmeas-vth)+(ub/tox/tox)*(xmeas-vth)**2)
    return sum( ((ymeas-y))**2 )
    
X =[400., 0.01, 0.01, 0.25]
X =[5., 1, 1, 5]

opts = {
      #'seed':1234, 
	'scaling_of_variables':[1e2, 1e-2, 1e-2, 1e-1], 
	'bounds':[[320, 0, 0, 0],[1000, 0.5, 0.5, 1.0]],
        'tolfun':1e-6,
	'tolx':1e-6,
	}


#print cma.fmin(f, X, 1, options=opts)
#exit()
for k, v in opts.iteritems():
    print k, "\t\t= ", v

es = cma.CMAEvolutionStrategy(X, 1, opts)
while not es.stop(): 
    X = es.ask()
    #for x in X:
    #    print x
    #raw_input()	
    es.tell(X, [f(x) for x in X])
    es.disp()
print('termination:', es.stop())
#cma.pprint(es.best.__dict__)
print es.best.__dict__['x']
print es.best.__dict__['f']





def f(x, u0, ua, ub, vth):
    tox = 1e-1
    return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
ymeas = array([f(x, 300, 0.1, 0.1, 0.2) for x in xmeas])

def f(x, u0, ua, ub, vth):
    tox = 1e-1
    return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)
yext0 = array([f(x, *array([500., 0.01, 0.01, 0.5])) for x in xmeas])

def f(x, u0, ua, ub, vth):
    tox = 1e-1
    return u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2)

   
X = es.best.__dict__['x'] 
yext = array([f(x, *X) for x in xmeas])
def f(X):
    [u0, ua, ub, vth] = X
    tox = 1e-1 
    y = u0/(1.0+(ua/tox)*(xmeas-vth)+(ub/tox/tox)*(xmeas-vth)**2)
    return sum( ((ymeas-y))**2 )



import matplotlib.pyplot as plt
plt.plot(xmeas, ymeas, 'r+')
plt.plot(xmeas, yext0, 'r')
plt.plot(xmeas, yext, 'g')
#plt.savefig('toto.png')
plt.show()
