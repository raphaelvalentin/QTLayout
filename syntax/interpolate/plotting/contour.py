#!/usr/local/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import SubplotParams
from matplotlib.ticker import MaxNLocator
import os
from function import *
from xplot import *

class contour(plot):
    def __init__(self, **kwargs):
        plot.__init__(self, **kwargs)

    def setProperties(self):
        plot.setProperties(self)
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['ytick.direction'] = 'out'
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	

    def plot(self, x, y, z):
        if not(len(z) == len(x) == len(y)):
            raise Exception('contour plot: x, y, z sizes are not equal.')
        size = (len(distinct(x)), len(distinct(y)))
        if not(size[0]*size[1] == len(z)):
            raise Exception('contour plot: x, y, z sizes are not coherent.')
	x = np.array(x); y = np.array(y); z = np.array(z)
	x.resize(size); y.resize(size); z.resize(size)
        self['plots'] = [(x, y, z, 10, {'color':'r'})]

    def getFigure(self):
        if not os.path.exists(self['filename']):
            self.setProperties()
            [(x, y, z, i, color)] = self['plots']
	    color['colors'] = 'r'
   	    CS = plt.contour(x, y, z, i, **color)
	    plt.clabel(CS, fontsize=19, inline=1, fmt='%g', colors='r')
            plt.draw()
            plt.savefig(self['filename'], dpi=self['dpi'])
        return self


	
