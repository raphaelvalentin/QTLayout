from math import log10
from libarray import *


def space(start, step1, stop, step2):
    step = step1/start+1
    start = log10(start)
    stop = log10(stop)
    step = log10(step)
    num = (stop-start)/step+1
    y = [10**(i*step+start) for i in xrange(int(num))]
    return y

    
if __name__=='__main__':
    x = 10.5
    
    print array(space(0+x, 0.5, 10+x, 1.0))-x


