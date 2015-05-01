# -*- coding: cp1252 -*-
from math import exp, sqrt
from libsolve import *
from libarray import *
import liboptim


txt = """-0.5	108.70
-0.45	120.00
-0.4	132.74
-0.35	147.06
-0.3	163.04
-0.25	180.72
-0.2	200.00
-0.15	220.59
-0.1	241.94
-0.05	263.16
0	283.02
0.05	300.00
0.1	312.50
0.15	319.15
0.2	319.15
0.25	312.50
0.3	300.00
0.35	283.02
0.4	263.16
0.45	241.94
0.5	220.59
0.55	200.00
0.6	180.72
0.65	163.04
0.7	147.06
0.75	132.74
0.8	120.00
0.85	108.70
0.9	98.68
0.95	89.82
1	81.97
"""

datax=array([])
datay=array([])
for i in xrange(0,len(txt.split()),2):
    datax.append(float(txt.split()[i]))
    datay.append(float(txt.split()[i+1]))

"""df * (dx) = -f"""
"""[df_da df_db] * [[da],[db]] = [-f]"""

##################
# wrap functions to levmar class

def f_(x, y0, X):
    u0, ua, ub = X[0], X[1], X[2]
    vth = X[3]
    tox = 1e-1
    return ((y0 - u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2))/y0)


def f(obj):
    global datax, datay
    global f0
    obj.function_call+=1
    X = array(obj.X.transpose()[0])
    f0 = f_(datax, datay, X)
    return array([f0.array]).transpose()

def J(obj):
    """ numerical jacobian matrix
    """
    global f0, f_, datax, datay, delta
    X = array(obj.X.transpose()[0])
    obj.function_call+=len(obj.X)
    delta=array([1e-6, 1e-6, 1e-6, 1e-6])

    row=array([])
    for i in xrange(len(X)):
        Z=X[:]
        Z[i]+=delta[i]
        row.append(-(f0 - f_(datax, datay, Z))/(delta[i]))
    J = row.transpose()

    del f0
    return J

#X =array([1000., 3e-9, 2e-18])
#X =array([4000., 0.5e-9, 0.5e-17])
X =array([400., 0.1, 0.04, 0.25])


#X =array([400., 0.1, 0.04, 0.25])

opt1 = liboptim.LMA(X, f, J)
opt1.echo = liboptim.echo
X1, norm_f, norm_dX = opt1.run()
#print X1




