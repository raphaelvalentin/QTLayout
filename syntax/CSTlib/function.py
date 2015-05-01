from math import *
from exceptions import Exception
import sys
from libarray import *
import re

# float in 64-bit coded
PRECISION_DECIMAL = 15
nan = float('NaN')

def sign(x):
    if x>=0: return 1
    return -1

def trunc(x, ndigits=PRECISION_DECIMAL):
    """ trunc(number/sequence, ndigits=PRECISION_DECIMAL) -> float/list point number(s) """
    def __trunc__(x):
        if x<>0:
            return round(x, -int(round(log10(abs(x)))+1)+ndigits)
	else:
	    return 0.0
    if hasattr(x, "__iter__"):
        return 	map(__trunc__, x)
    else: 
        return __trunc__(x)


def linspace(start, stop, num=0, step=None):
    start = float(start)
    stop = float(stop)
    if step == None:
        if num == 1:
	     step = stop-start
	else:
    	     step = (stop-start)/(float(num)-1.0)
    else:
        step = float(step)
	if step == 0.0:
	    num = 0
	else:
	    num = (stop-start)/step+1.0
    return array(trunc([i*step+start for i in xrange(int(num + pow(2,-45)))]))


def logspace(start, stop, num=1, step=None, prec=10):
    if start<=0.0 or stop<=0.0:
        return []
    if step == None:
        y = linspace(log10(start), log10(stop), num=num)
    else:
        if step<=0:
	    return []
        y = linspace(log10(start), log10(stop), step=log10(step))
    return array(trunc([10**i for i in y], prec))


def derivative(x, y, n=1):
    def __derivative__(x, y):
        xr = []
        yr = []
        for i in xrange(len(x)-1):
            xr.append(((x[i]+x[i+1])*0.5))
            yr.append(((y[i]-y[i+1])/(x[i]-x[i+1])))
        return (xr,yr)
    if len(x)<n+1:
        raise Exception('Warning : not enough points for applying the derivative.')
    for i in xrange(n):
        (x, y) = __derivative__(x, y)
    return (x, y)


def flatten(sequence):
    """ yield each element of an irregular list (or tuple, dict...->__instance)"""
    for el in sequence:
        if isinstance(el, (list, tuple)):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def rect(r, phi):
    return r*(cos(phi)+1j*sin(phi))


def closest1(flist, value):
    if len(flist) > 0:
        b1 = flist[0]
        for v in flist:
            if abs(v-value)<abs(b1-value):
                b1 = v
        return b1
    else:
        raise Exception('Warning : list length is too short')

nearest1 = closest1

def closest2(flist, value):
    if len(flist) > 1:
        b1 = flist[0]
        b2 = flist[1]
        for v in flist:
            if abs(v-value)<abs(b1-value) and v<>b2:
                b1 = v
        for v in flist:
            if abs(v-value)<abs(b2-value) and b1<>v:
                b2 = v
        return b1, b2
    else:
        raise Exception('Warning : list length is too short')

nearest2 = closest2

def interpolate(x, y, x0):
    xp0, xp1 = closest2(x, x0)
    i0, i1 = x.index(xp0), x.index(xp1)
    a = (y[i1]-y[i0])/(x[i1]-x[i0])
    b = y[i0]-a*x[i0]
    y0 = (a*x0+b)
    return y0
    
    
def transpose(m):
    return zip(*m)
    
    
def distinct(seq):
    result = []
    [result.append(i) for i in seq if not result.count(i)]
    return result


def meshgrid(x, y):
    _x = array([])
    _y = array([])
    for j in y:
        for i in x:
            _x.append(i)
            _y.append(j)
    return (_x, _y)

def unit(s):
    """Takes a string and returns the equivalent float.
    '3.0u' -> 3.0e-6"""
    mult={'t':1.0e12,
          'g':1.0e9,
          'meg':1.0e6,
          'k':1.0e3,
          'mil':25.4e-6,
          'm':1.0e-3,
          'u':1.0e-6,
          'n':1.0e-9,
          'p':1.0e-12,
          'f':1.0e-15,
	  'a':1.0e-18}
    m=re.search('^([0-9e\+\-\.]+)(t|g|meg|k|mil|m|u|n|p|f)?',s.lower())
    if m.group(2):
        return float(m.group(1))*mult[m.group(2)]
    else:
        return float(m.group(1))



def gradspace(start, step_beg, stop, step_end, maxnodes=1e3):
    a = (step_end-step_beg)/(stop-start)
    b = step_beg-a*start
    r = array([start])
    x = start
    for i in xrange(int(maxnodes)):
        step = x * a + b
        x = x + step
	if x>=stop:
	    break
        r.append(x)
    if i == maxnodes-1:
        raise Exception('max nodes reach')
    if (stop-r[-1])/step_end<0.2:
        r.append(stop)
        return r[:-1]
    r.append(stop)
    return r

_gradspace = gradspace

def gradspace2(start, step_beg, stop, step_end):
    r1 = _gradspace(start, step_beg, stop, step_end)
    x1, y1 = step_end, r1[-1]
    r = _gradspace(start, step_beg, stop, step_end*0.99)
    x2, y2 = step_end*0.99, r[-1]
    a = (y2-y1)/(x2-x1)
    b = y1-a*x1
    x = (stop-b)/a
    r = _gradspace(start, step_beg, stop, x)
    assert r[-1]-r[-2]>0.9*x
    r[-1] = stop
    return r


def gradspace1(start, step, stop, maxstep=None,  minstep=None, expansion=1.3):
    if not maxstep:
        maxstep=1e100
    if not minstep:
        minstep=0.0
    step = (step/abs(stop-start))**(1./expansion)
    y = [start]
    i = 1
    while i<1e3:
        s = ((i*step)**expansion) * (stop-start) + start
	if expansion>1.0 and abs(s-y[-1])>maxstep:
	    s = y[-1] + sign(stop-start)*maxstep
	elif expansion<1.0 and abs(s-y[-1])<minstep:
	    s = y[-1] + sign(stop-start)*minstep
	if (stop>start and s>stop) or (stop<start and s<stop):
	    break
        y.append(s)
	i = i+1
    y.append(stop)	
    return y

def integ(f,X0,X1,n):
#integration numerique a partir d une fonction
    H=X1-X0
    p=H/n
    I=0.
    h=p/4.
    x0=X0
    while (x0+p)<X1:
        x1=x0+h
        x2=x1+h
        x3=x2+h
        x4=x3+h
        I+=14.*f(x0)+64.*f(x1)+24.*f(x2)+64.*f(x3)+14.*f(x4)
        x0=x4
    return I*h/45.
