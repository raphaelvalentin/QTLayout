from xyplot3 import *
from libarray import *

a = array([1,2,3,4,5,6,7,8,9,10])
b = sin(a/4.)

plt = plot()
plt.plot(a, b)
plt.savefig(filename='toto.png', force=True)
