from CSTlib import *
from math import *

def Octogon(diameter=1.0):
    octogon1 = Primitive()
    r = 0.5*diameter
    a = diameter/(1.0+sqrt(2.0))
    octogon1.append( Point(-0.5*a, -r, 0) )
    octogon1.append( Point(0.5*a, -r, 0) )
    octogon1.append( Point(r, -0.5*a, 0) )
    octogon1.append( Point(r, 0.5*a, 0) )
    octogon1.append( Point(0.5*a, r, 0) )
    octogon1.append( Point(-0.5*a, r, 0) )
    octogon1.append( Point(-r, 0.5*a, 0) )
    octogon1.append( Point(-r, -0.5*a, 0) )
    return octogon1

octogon = Octogon
