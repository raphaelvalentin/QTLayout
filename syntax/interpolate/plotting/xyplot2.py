#!/usr/local/bin/python

import matplotlib.pyplot as plt
from function import *
from libarray import *

x = linspace(1, 10, 50)
y = sin(x)

plt.figure(1)
plt.subplot(211)
plt.plot(x, y, 'r', label='toto')
plt.plot(x, y, 'b', label='toto1')
plt.legend()
plt.subplot(212)

plt.plot(x,y, 'b', label='tat')
plt.legend()
plt.show()
