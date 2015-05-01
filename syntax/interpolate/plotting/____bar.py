#!/usr/local/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import SubplotParams
from matplotlib.ticker import MaxNLocator
import os
from xplot import *

class bar(plot):

    def __init__(self, **kwargs):
        plot.__init__(self, **kwargs)
        self.xlogscale = kwargs.pop('xlogscale', False)
        self.ylogscale = kwargs.pop('ylogscale', False)
        self.width = kwargs.pop('width', 0.25)

    def plot(self, *args, **kwargs):
	for xy in args:
            x, y = xy
	    if not isinstance(x, (list, tuple)):
	        raise Exception('xyplot : x is not a list or a tuple.')
	    if not isinstance(y, (list, tuple)):
	        raise Exception('xyplot : y is not a list or a tuple.')
	    if not len(x) == len(y):
	        raise Exception('xyplot : x and y sizes are different.')
            if 'label' in kwargs and True in map(lambda elt: elt[2].get('label', '')==kwargs['label'], self['plots']):
                del kwargs['label']
            properties = {'linewidth':2 }
            properties.update(**kwargs)
            self['plots'].append((x, y, properties))

    def setProperties(self):
        plot.setProperties(self)

    def getFigure(self):
        if not os.path.exists(self['filename']):
            self.setProperties()
            legend = False
	    for (x, y, parameters) in self['plots']:
    	        plt.bar(x, y, **parameters)
                if 'label' in parameters: legend = True
            if legend:
                plt.legend(loc=0, prop={'size':self.fontsize})
     	        # transparent legend
                leg = self.ax.legend(loc='best', fancybox=False)
                leg.get_frame().set_alpha(0.5)
            plt.draw()
            plt.savefig(self['filename'], dpi=self['dpi'])
	    
        return self
	
