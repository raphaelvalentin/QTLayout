import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import SubplotParams
from matplotlib.ticker import MaxNLocator
import os

class plot(dict):
    def __init__(self, **kwargs):
        self['filename'] = kwargs.get('filename', 'image.png')
        self['dpi'] = kwargs.get('dpi', 100)
        self['caption'] = kwargs.get('caption', '')
        self['figsize'] = kwargs.get('figsize', (8, 6))
        self['plots'] = kwargs.get('plots', [])
        self.xlabel = kwargs.get('xlabel', '')
        self.ylabel = kwargs.get('ylabel', '')
        self.fontsize = kwargs.get('fontsize', 19)
	self.xlogscale = False
	self.ylogscale = False

    def setProperties(self):
	plt.clf()
        plt.rc('font', family='sans-serif', weight='big', size=self.fontsize)
        plt.rc('figure', figsize=(8,6))
        plt.rc('figure', dpi=self['dpi'])
        plt.rc('figure.subplot', left=0.20, bottom=0.13, right=0.93, top=0.93, wspace=0.001, hspace=0.1)
        plt.rc('lines', markersize=6)
        plt.rc('axes', labelsize=self.fontsize)
        plt.rc('axes', color_cycle=('red', 'blue', 'green', 'black', 'grey', 'yellow'))
        plt.rc('axes', grid=True)
        plt.rc('xtick.major', size=8)           # major tick size in points
        plt.rc('xtick.minor', size=5)           # minor tick size in points
        plt.rc('xtick.major', width=1.5)        # major tick width in points
        plt.rc('xtick.minor', width=1.5)        # minor tick width in points
        plt.rc('xtick.major', pad=4)            # distance to major tick label in points
        plt.rc('xtick.minor', pad=4)            # distance to the minor tick label in points
        plt.rc('xtick', color='k')                # color of the tick labels
        plt.rc('xtick', labelsize=self.fontsize)       # fontsize of the tick labels
        plt.rc('xtick', direction='in')           # direction: in, out, or inout
        plt.rc('ytick.major', size=8)           # major tick size in points
        plt.rc('ytick.minor', size=5)           # minor tick size in points
        plt.rc('ytick.major', width=1.5)        # major tick width in points
        plt.rc('ytick.minor', width=1.5)        # minor tick width in points
        plt.rc('ytick.major', pad=4)            # distance to major tick label in points
        plt.rc('ytick.minor', pad=4)            # distance to the minor tick label in points
        plt.rc('ytick', color='k')                # color of the tick labels
        plt.rc('ytick', labelsize=self.fontsize)       # fontsize of the tick labels
        plt.rc('ytick', direction='in')           # direction: in, out, or inout

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)

	# max number of ticks
	plt.locator_params(nbins=7)
	
        if self.xlogscale:
            self.ax.set_xscale('log')
        if self.ylogscale:
            self.ax.set_yscale('log')
        ax_r = plt.gca() #for each axis or whichever axis you want you should

        #plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))

from xyplot31 import plot

__all__ = ['plot']    
