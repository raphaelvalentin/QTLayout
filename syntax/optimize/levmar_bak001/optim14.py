# -*- coding: cp1252 -*-
from math import exp, sqrt
from libsolve import *
from libarray import *
import lm


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

def f(X):
    global datax, datay
    def f_(x, y0, X):
        u0, ua, ub = X[0], X[1], X[2]
        vth = X[3]
        tox = 1e-1
        return ((y0 - u0/(1.0+(ua/tox)*(x-vth)+(ub/tox/tox)*(x-vth)**2))/y0)
    return f_(datax, datay, X)


X =array([320., 0.1, 0.04, 0.25])
opt1 = lm.LMA(f, X, xtol=1e-3)
X, norm_f, norm_dX = opt1.run()




