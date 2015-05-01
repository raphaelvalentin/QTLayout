from xplot import *
from function import *
from libarray import *

a = xyplot(xlogscale=True)
x = linspace(0,10, step=0.1)
y = sin(x)

a.plot((x,y), color='b', ls='.')

print a.getFigure()['filename']
