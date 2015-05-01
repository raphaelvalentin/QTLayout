from libarray import *

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
	    num = (stop-start)/step+1
    return array([trunc(i*step+start) for i in xrange(int(num))])


def smoothspace(start, step_beg, stop, step_end):
    a = (step_end-step_beg)/((stop-step_end)-start)
    b = step_beg-a*start
    r = array([start])
    x = start
    while x<stop:
        step = x * a + b
	x = x + step
        r.append(x)
    return r
    
for i in smoothspace(0.0, 0.1, 2.0, 0.2)    :
    print i, 0.0
    print i, 1.0
    print 
    

    
