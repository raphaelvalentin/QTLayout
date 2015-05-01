from libarray import shape
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import SubplotParams
from matplotlib.ticker import MaxNLocator
from os.path import isfile

class plot(object):
    def __init__(self, **kwargs):
        self.dpi = kwargs.get('dpi', 100)
        self.caption = kwargs.get('caption', '')
        self.figsize = kwargs.get('figsize', (8, 6))
        self.xlabel = kwargs.get('xlabel', '')
        self.ylabel = kwargs.get('ylabel', '')
        self.fontsize = kwargs.get('fontsize', 19)
        self.xlogscale = kwargs.pop('xlogscale', False)
        self.ylogscale = kwargs.pop('ylogscale', False)
	self._labels = list()
	self._plots = list()
        self.filename = kwargs.get('filename', 'image.png')
	
	#plt.clf()
        plt.rc('font', family='sans-serif', weight='big', size=self.fontsize)
        plt.rc('figure', figsize=(8,6))
        plt.rc('figure', dpi=self.dpi)
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

    def plot(self, *args, **kwargs):
        if not 'linewidth' in kwargs:
	    kwargs['linewidth'] = 2
        i = 0
        while i<len(args):
	    x = args[i]
	    i = i+1
	    if not isinstance(x, list):
	        raise Exception('x is not a list')
	    if len(shape(x)) == 1:
	        y = args[i]
		i = i+1
	    elif len(shape(x)) == 2:
	        x, y = x
	    elif len(shape(x)) == 3:
	        raise Exception('the shape of x is not correct')
   	    if not isinstance(y, list):
	        raise Exception('y is not a list')
	    if len(shape(y))<>1:
	        raise Exception('the shape of y is not 1')
	    if len(x)<>len(y):
	        raise Exception('x and y sizes are different.')
	    if 'label' in kwargs:
		if kwargs['label'] in self._labels:
		    del kwargs['label']
		else:
		    self._labels.append(kwargs['label'])
	    if i<len(args) and isinstance(args[i], str):
		linespec = args[i]
		i = i+1
		self._plots.append(plt.plot(x, y, linespec, **kwargs))
	    else:
		self._plots.append(plt.plot(x, y, **kwargs))

    def scatter(self, *args, **kwargs):
        markeredgecolor = kwargs.pop('color', 'r')
	markersize = kwargs.pop('size', 6)
        properties = {'marker':'s', 'markersize':markersize, 'linewidth':0, 
                      'markerfacecolor':'none', 'markeredgecolor':markeredgecolor, 'markeredgewidth':2
                     }
	properties.update(**kwargs)
	self.plot(*args, **properties)
        
    def show(self):
        if len(self._labels):
            plt.legend(loc=0, prop={'size':self.fontsize})
     	    # transparent legend
            leg = self.ax.legend(loc='best', fancybox=False)
            leg.get_frame().set_alpha(0.5)  
        plt.draw()
	plt.show()

    def getFigure(self, force=False):
        if len(self._plots)>0 and (force or not(isfile(self.filename))):
            if len(self._labels):
                plt.legend(loc=0, prop={'size':self.fontsize})
     	        # transparent legend
                leg = self.ax.legend(loc='best', fancybox=False)
                leg.get_frame().set_alpha(0.5)  
            plt.draw()
	    plt.savefig(self.filename, dpi=self.dpi)
	return {'filename':self.filename, 'caption':self.caption}

    def savefig(self, filename, dpi=100, force=False):
        self.filename = filename
	self.dpi = dpi
        self.getFigure(force=force)

xyplot = plot
        

if __name__ == '__main__':

    from libarray import *
    a = array([1,2,3,4,5,6,7,8,9,10])
    b = sin(a)
    c = a+0.5
    d = cos(c)

    x = plot(xlabel='x', ylabel='y')
    x.plot([a, b], color='r', label='toto')
    x.scatter([c, d], color='r', label='toto')
    x.savefig('totox.png')
    

		
	    


    
