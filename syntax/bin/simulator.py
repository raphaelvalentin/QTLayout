from subprocess import Popen, PIPE, STDOUT

class fasthenry(object):
    def __init__(self, netlist=''):
        self.netlist = netlist

    def run(self, verbose=False, refinement=1):
        filename = 'source.inp'
        with open(filename, 'w') as f:
            f.write(str(self.netlist))
        cmd = ['fasthenry.exe', filename]
        if refinement>1:
            cmd = ['fasthenry.exe', '-i', str(refinement), filename]
        proc = Popen(" ".join(cmd), stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True)
        stdout = proc.communicate()[0]
        if verbose:
            print stdout
        self.raw = self.parse(stdout)

    def parse(self, stdout):
        lines = stdout.split('\n')
        raw = {'freq':[]}
        isIMPEDANCE_MATRIX = False
        indxLine = 0
        while indxLine<len(lines):
            line = lines[indxLine]
            indxLine += 1
            if 'Frequency' in line:
                frequency = float(line.strip().split()[-1])
                raw['freq'].append(frequency)
                continue
            if 'Impedance matrix' in line:
                isIMPEDANCE_MATRIX = True
                s = line.split()
                size = int(s[s.index('x')-1]), int(s[s.index('x')+1])
                if size == (1,1):
                    if not raw.has_key('z11'):
                        raw['z11'] = []
                if size == (2,2):
                    if not raw.has_key('z11'):
                        raw['z11'] = []
                    if not raw.has_key('z12'):
                        raw['z12'] = []
                    if not raw.has_key('z21'):
                        raw['z21'] = []
                    if not raw.has_key('z22'):
                        raw['z22'] = []
                continue
            if isIMPEDANCE_MATRIX:
                if size == (1,1):
                    s = line.replace('j', '').split()
                    d = [float(x) for x in s]
                    raw['z11'].append(complex(d[0], d[1]))
                if size == (2,2):
                    s = line.replace('j', '').split()
                    d = [float(x) for x in s]
                    raw['z11'].append(complex(d[0], d[1]))
                    raw['z12'].append(complex(d[2], d[3]))
                    line = lines[indxLine]
                    indxLine += 1
                    s = line.replace('j', '').split()
                    d = [float(x) for x in s]
                    raw['z21'].append(complex(d[0], d[1]))
                    raw['z22'].append(complex(d[2], d[3]))
                isIMPEDANCE_MATRIX = False
                continue
        return raw                        
            

class FastHenry(list):
    def __init__(self):
        pass
    def __str__(self):
        self.append('.end')
        return "\n".join([str(line) for line in self])

class Title:
    def __init__(self, title='My Fast Henry Script'):
        self.title = title
    def __str__(self):
        return "** %r" % self.title
    
class Units:
    units = 'um'
    def __init__(self, units='um'):
        Units.units = units
    def __str__(self):
        return ".Units %r" % Units.units

class Default:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
    def __str__(self):
        s = ['.Default']
        if self.kwargs.has_key('z'):
            s.append( "z=%r"%self.kwargs['z'] )
        if self.kwargs.has_key('sigma'):
            if Units.units=='um':
                s.append( "sigma=%r"%(self.kwargs['sigma']*1e-6) )
            if Units.units=='mm':
                s.append( "sigma=%r"%(self.kwargs['sigma']*1e-3) )
        return " ".join(s)
        
class Node:
    def __init__(self, x, y, z=None):
        self.x = x
        self.y = y
        self.z = z
        self.name = 'N%d' % Node.indx
        Node.indx += 1
    def __str__(self):
        if self.z:
            return "%s x=%r y=%r z=%r" % (self.name, self.x, self.y, self.z)
        else:
            return "%s x=%r y=%r" % (self.name, self.x, self.y)
Node.indx = 1
        
class Segment:
    def __init__(self, node1, node2, w=1, h=1, **kwargs):
        self.node1 = node1
        self.node2 = node2
        self.name = 'E%d' % Segment.indx
        Segment.indx += 1
        self.kwargs = kwargs
        self.w = w
        self.h = h
    def __str__(self):
        s = [self.name, self.node1.name, self.node2.name, 'w=%r'%self.w, 'h=%r'%self.h]
        if len(self.kwargs):
            s.extend( ["%s=%r"%(key, value) for key, value in self.kwargs.iteritems()] )
        return " ".join(s)
Segment.indx = 1
    
class Port:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
    def __str__(self):
        s = ['.external', self.node1.name, self.node2.name]
        return " ".join(s)

class Freq:
    def __init__(self, fmin, fmax, ndec=1):
        self.fmin = fmin
        self.fmax = fmax
        self.ndec = ndec
    def __str__(self):
        s = ['.freq', 'fmin=%r'%self.fmin, 'fmax=%r'%self.fmax, 'ndec=%r'%self.ndec]
        return " ".join(s)
    
        
        
if __name__=='__main__':
    netlist = FastHenry()
    netlist.append( Title('my script') )
    netlist.append( Default(z=0, sigma=5e7) )
    netlist.append( Freq(1e6, 1e9, 1) )
    n1 = Node(x=0, y=0)
    n2 = Node(x=1000, y=0)
    n3 = Node(x=1000, y=1000)
    n4 = Node(x=0, y=1000)
    n5= Node(x=0, y=1)
    n6 = Node(x=200, y=200)
    n7 = Node(x=800, y=200)
    n8 = Node(x=800, y=800)
    n9 = Node(x=200, y=800)
    n10 = Node(x=200, y=201)

    netlist.extend([n1,n2,n3,n4,n5,n6,n7,n8,n9,n10])

    width = 10
    thickness = 3.4
    netlist.append( Segment(n1,n2,w=width,h=thickness, nhinc=3, nwinc=3) )
    netlist.append( Segment(n2,n3,w=width,h=thickness, nhinc=3, nwinc=3) )
    netlist.append( Segment(n3,n4,w=width,h=thickness, nhinc=3, nwinc=3) )
    netlist.append( Segment(n4,n5,w=width,h=thickness, nhinc=3, nwinc=3) )
    
    width = 8
    thickness = 3.4
    netlist.append( Segment(n6,n7,w=width,h=thickness, nhinc=3, nwinc=3) )
    netlist.append( Segment(n7,n8,w=width,h=thickness, nhinc=3, nwinc=3) )
    netlist.append( Segment(n8,n9,w=width,h=thickness, nhinc=3, nwinc=3) )
    netlist.append( Segment(n9,n10,w=width,h=thickness, nhinc=3, nwinc=3) )

    netlist.append( Port(n1, n5) )
    netlist.append( Port(n6, n10) )

    
    simu = fasthenry(netlist)
    simu.run(verbose=False)
    print simu.raw
