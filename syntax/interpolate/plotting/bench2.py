from xplot import *
from function import *
from libarray import *

a = contour()
x = linspace(0,10, step=0.1)
y = linspace(0,10, step=0.1)

X, Y = meshgrid(x, y)

Z = sin(X)*sin(Y)

a.plot(X, Y, Z)

print a.getFigure()['filename']
