from libarray import shape
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import SubplotParams
from matplotlib.ticker import MaxNLocator
from os.path import isfile

class plot(object):
    __properties__ = {'figsize':(8, 6),
                      'xlabel':'',
                      'ylabel':'',
                      'xlogscale':False,
                      'ylogscale':False,
                      'dpi':100,
                      'filename':'image.png',
                      'force':False,
                     }

    def __init__(self, **properties):
        self.properties = dict(plot.__properties__)
        self.properties.update(properties)
        self.lines = list()

    def xlabel(self, xlabel):
        self.properties['xlabel'] = xlabel
    def ylabel(self, ylabel):
        self.properties['ylabel'] = ylabel
    def force(self, force):
        self.properties['force'] = force
	
    def plot(self, *args, **properties):
        if not 'linewidth' in properties:
	    properties['linewidth'] = 2
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
	    if i<len(args) and isinstance(args[i], str):
		linespec = args[i]
		i = i+1
		self.lines.append((x, y, linespec, properties))
	    else:
		self.lines.append((x, y, properties))
		
    def scatter(self, *args, **properties):
        markeredgecolor = properties.pop('color', 'r')
	markersize = properties.pop('size', 6)
        _properties = {'marker':'s', 'markersize':markersize, 'linewidth':0, 
                      'markerfacecolor':'none', 'markeredgecolor':markeredgecolor, 'markeredgewidth':2
                      }
	_properties.update(**properties)
	self.plot(*args, **_properties)
		
    def savefig(self, filename, dpi=100, force=False):
        self.properties['filename'] = str(filename)
	self.properties['dpi'] = int(dpi)
	self.properties['force'] = bool(force)
	
	if len(self.lines)==0:
            raise Exception('no line has been defined in the plot.')
        if isfile(self.properties['filename']) and self.properties['force']==False:
            print 'Warning: The image file \'%s\' has not been updated.'%self.properties['filename']
            return
	
        plt.rc('font', family='sans-serif', weight='big', size=19)
        plt.rc('figure', figsize=self.properties['figsize'])
        plt.rc('figure', dpi=self.properties['dpi'])
        plt.rc('figure.subplot', left=0.20, bottom=0.13, right=0.93, top=0.93, wspace=0.001, hspace=0.1)
        plt.rc('lines', markersize=6)
        plt.rc('axes', labelsize=19)
        plt.rc('axes', color_cycle=('red', 'blue', 'green', 'black', 'grey', 'yellow'))
        plt.rc('axes', grid=True)
        plt.rc('xtick.major', size=8)           # major tick size in points
        plt.rc('xtick.minor', size=5)           # minor tick size in points
        plt.rc('xtick.major', width=1.5)        # major tick width in points
        plt.rc('xtick.minor', width=1.5)        # minor tick width in points
        plt.rc('xtick.major', pad=4)            # distance to major tick label in points
        plt.rc('xtick.minor', pad=4)            # distance to the minor tick label in points
        plt.rc('xtick', color='k')                # color of the tick labels
        plt.rc('xtick', labelsize=19)       # fontsize of the tick labels
        plt.rc('xtick', direction='in')           # direction: in, out, or inout
        plt.rc('ytick.major', size=8)           # major tick size in points
        plt.rc('ytick.minor', size=5)           # minor tick size in points
        plt.rc('ytick.major', width=1.5)        # major tick width in points
        plt.rc('ytick.minor', width=1.5)        # minor tick width in points
        plt.rc('ytick.major', pad=4)            # distance to major tick label in points
        plt.rc('ytick.minor', pad=4)            # distance to the minor tick label in points
        plt.rc('ytick', color='k')                # color of the tick labels
        plt.rc('ytick', labelsize=19)       # fontsize of the tick labels
        plt.rc('ytick', direction='in')           # direction: in, out, or inout

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        plt.xlabel(self.properties['xlabel'])
        plt.ylabel(self.properties['ylabel'])

	# max number of ticks
	plt.locator_params(nbins=7)
	
        if self.properties['xlogscale']:
            self.ax.set_xscale('log')
        if self.properties['ylogscale']:
            self.ax.set_yscale('log')
        ax_r = plt.gca() #for each axis or whichever axis you want you should
        
        labels = []
        for line in self.lines:
            if 'label' in line[-1]:
               if line[-1]['label'] in labels:
                   del line[-1]['label'] 
               else: 
                   labels.append(line[-1]['label'])
                   
        for line in self.lines:
            if len(line)==4:
                x, y, linespec, properties = line
                plt.plot(x, y, linespec, **properties)
            elif len(line)==3:
                x, y, properties = line
                plt.plot(x, y, **properties)
                 
        if len(labels):
           plt.legend(loc=0, prop={'size':19})
           # transparent legend
           leg = self.ax.legend(loc='best', fancybox=False)
           leg.get_frame().set_alpha(0.5)  
        plt.draw()
        plt.savefig(self.properties['filename'], dpi=self.properties['dpi'])



if __name__ == '__main__':
    from libarray import *
    from function import *

    x = linspace(0, 10, step=0.2)
    y1 = sin(x)
    y2 = cos(x)
    y3 = atan(x)

    plt1 = plot(xlabel='x', ylabel='y1')
    plt2 = plot(xlabel='x', ylabel='y2, y3')

    plt1.plot(x, y1, color='r', label='sin(x)')
    plt2.scatter(x, y2, color='r', label='sin(x)')
    plt2.plot(x, y3, color='b', label='atan(x)')

    plt1.savefig(filename='image1.png', force=True)
    plt2.savefig(filename='image2.png', force=True)

		


