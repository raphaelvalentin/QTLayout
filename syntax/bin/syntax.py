from math import pi, sqrt

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
        return ".Units %s" % Units.units

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
    def __init__(self, node1, node2, w=1, h=1, nwinc=1, nhinc=1, **kwargs):
        self.node1 = node1
        self.node2 = node2
        self.name = 'E%d' % Segment.indx
        Segment.indx += 1
        self.kwargs = kwargs
        self.w = w
        self.h = h
        self.nwinc=nwinc
        self.nhinc=nhinc
    def __str__(self):
        s = [self.name, self.node1.name, self.node2.name, 'w=%r'%self.w, 'h=%r'%self.h, 'nhinc=%r'%self.nhinc, 'nwinc=%r'%self.nwinc ]
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

class Equiv:
    def __init__(self, *nodes):
        self.nodes = nodes
    def __str__(self):
        s = ['.equiv'] + [node.name for node in self.nodes]
        return " ".join(s)

class Freq:
    def __init__(self, fmin, fmax, ndec=1):
        self.fmin = fmin
        self.fmax = fmax
        self.ndec = ndec
    def __str__(self):
        s = ['.freq', 'fmin=%r'%self.fmin, 'fmax=%r'%self.fmax, 'ndec=%r'%self.ndec]
        return " ".join(s)
    
        
        
