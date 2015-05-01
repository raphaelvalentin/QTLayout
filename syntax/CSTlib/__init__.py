__all__ = ['Brick', 'Copy', 'Component', 'CST', 'Layers','Layer',
           'DiscretePort', 'Boundary',
           'Unit', 'Material', 'Mesh', 'Solver', 'Extrude',
           'Solid', 'Solids', 'Group', 'MeshSettings',
           'Transform', 'Subtract', 'Intersect', 'LumpedElement', 'Blend', 'walk2intersect'
           ]

from math import pi, log10
from formal import *
from exceptions import *

def isPoint(obj):
    try:
        obj.x
        obj.y
        return True
    except:
        return False


Angle = lambda pt: atan2(pt.y, pt.x)/pi*180.

def walk2intersect( shape1, shape2, start=0, end=0, increment=1):
    curve1 = Primitive()
    for i in xrange(start, end, increment):
        curve1.append( shape1[i] )
        for j in xrange(len(shape2)-1):
            segment1 = Segment(shape1[i], shape1[i+increment])
            segment2 = Segment(shape2[j], shape2[j+1])
            point1 = segment1.Intercept( segment2 )
            if point1:
                curve1.append( point1 )
                return curve1








from libarray import *
from linalg import solve
from function import linspace
from newton import NewtonRaphson as fmin
from math import atan2
Angle = lambda pt: atan2(pt.y, pt.x)/pi*180.

               
               
def Blend(curve, n=100, radius=1.0):
    assert isinstance(curve, list), 'TypeError: curve it not a Primitive'
    assert len(curve)==3, 'TypeError: number of point is different of 3'
    Angle = lambda pt: atan2(pt.y, pt.x)/pi*180.
    n = int(n)
    point1, point2, point3 = curve
    arr = lambda x: round(x*1e6)/1e6
    point1 = Point(arr(point1.x), arr(point1.y), arr(point1.z))
    point2 = Point(arr(point2.x), arr(point2.y), arr(point2.z))
    point3 = Point(arr(point3.x), arr(point3.y), arr(point3.z))
    
    isreversed = False
    if not(point1.x<=point2.x<=point3.x):
        point1, point3 = point3, point1
        isreversed = True
    line1 = Line(point1, point2)
    line2 = Line(point2, point3)
    if isinf(line1.slope): line1 = Line(point1, point2+Point(1e-6, 0, 0))
    if isinf(line2.slope): line2 = Line(point2, point3+Point(1e-6, 0, 0))
    def f(x):
        if point1.x<=point2.x<=point3.x:
            point4 = point2.Translate2(line=line1, radius=-x)
            point5 = point2.Translate2(line=line2, radius=x)
        line4 = line1.Rotate(center=point4, angle=(0,0,90))
        line5 = line2.Rotate(center=point5, angle=(0,0,90))
        point6 = line4.Intercept(line5)
        return (point4.x-point6.x)**2+(point4.y-point6.y)**2 - radius**2
    x = fmin(f=f, x0=radius, verbose=False).run()
    if point1.x<=point2.x<=point3.x:
        point4 = point2.Translate2(line=line1, radius=-x)
        point5 = point2.Translate2(line=line2, radius=x)
    line4 = line1.Rotate(center=point4, angle=(0,0,90))
    line5 = line2.Rotate(center=point5, angle=(0,0,90))
    point6 = line4.Intercept(line5)
    ang1 = Angle(point4-point6)
    ang2 = Angle(point5-point6)
    curve3 = Primitive()
    curve3.append( point1 )
    for theta in linspace(ang1, ang2, n):
        curve3.append( Point(radius*cos(theta*pi/180), radius*sin(theta*pi/180), 0)+point6 )
    curve3.append( point3 )
    if isreversed:
        return reversed(curve3)
    return curve3




def norm(M):
    err = 0.
    for i in M:
        for j in i:
            err+=j*j
    return math.sqrt(err)


class Brick(object):
    _names = {'Brick':1}
    def __init__(self, **kwargs):
        name = kwargs.get('name', 'Brick')
        if name in Brick._names:
            self.name = "{name}{i}".format(name=name, i=Brick._names[name])
            Brick._names[name] += 1
        else:
            self.name = "{name}".format(name=name)
            Brick._names[name] = 1
        self.component = kwargs.get('component', 'component1')
        self.material = kwargs.get('material', 'PEC')
        if not 'points' in kwargs:
            self.xrange = kwargs.get('xrange', (0, 1))
            self.yrange = kwargs.get('yrange', (0, 1))
            self.zrange = kwargs.get('zrange', (0, 1))
        else:
            self.points = kwargs.get('points')
            self.xrange = ( self.points[0].x, self.points[1].x)
            self.yrange = ( self.points[0].y, self.points[1].y)
            self.zrange = ( self.points[0].z, self.points[1].z)
    def __str__(self):
        return "\n".join([ 'With Brick',
                           '  .Reset',
                           '  .Name "{name}"'.format(name=self.name),
                           '  .Component "{component}"'.format(component=self.component),
                           '  .Material "{material}"'.format(material=self.material),
                           '  .Xrange "{xmin}", "{xmax}"'.format(xmin=self.xrange[0], xmax=self.xrange[1]),
                           '  .Yrange "{ymin}", "{ymax}"'.format(ymin=self.yrange[0], ymax=self.yrange[1]),
                           '  .Zrange "{zmin}", "{zmax}"'.format(zmin=self.zrange[0], zmax=self.zrange[1]),
                           '  .Create',
                           'End With',
                           ''
                        ])
    def __copy__(self):
        brick = Brick(component=self.component, material=self.material, xrange=self.xrange, yrange=self.yrange, zrange=self.zrange)
        return brick
        

def Copy(obj):
    if '__copy__' in dir(obj):
        return obj.__copy__()
    d = {}
    for k, v in obj.__dict__.iteritems():
        d[k] = v
    l = []
    if isinstance(obj, list):
        for o in obj:
            l.append(Copy(o))
    new = obj.__class__(*l)
    if 'name' in new.__dict__:
        d['name'] = new.__dict__['name']
    new.__dict__.update(d)
    return new
    

class Component(object):
    @staticmethod
    def New(name):
        return "\n".join(['Component.New "{name}"'.format(name=name), '' ])


class CST(list):
    Formatting = ":.5f"
    def __init__(self, **kwargs):
        self.component = kwargs.get('component', 'component1')
        self.name = kwargs.get('name', 'object1')
    def __str__(self):
        s = ['Sub Main ()', '']
        for line in self:
            if not isinstance(line, CST):
                s.append( str(line) )
            else:
                s.append( "\n".join(str(line).split('\n')[2:-1]) )
        s.extend(['End Sub'])
        return "\n".join(s)
            

class Layer(list):
    unit = 'A'
    prec = 6
    def __init__(self, name, material='', zmax=None, thickness=None, zmin=None, **kwargs):
        self.name = name
        self.material = material
        for k, v in kwargs.iteritems():
            self.__setattr__(k, v)
        if Layer.unit=='A':
            factor = 1e-4
        elif Layer.unit=='nm':
            factor = 1e-3
        elif Layer.unit=='um':
            factor = 1.0
        if zmax<>None and thickness:
            zrange = (zmax-thickness)*factor, zmax*factor
        elif zmin<>None and thickness:
            zrange = zmin*factor, (zmin+thickness)*factor
        elif zmin<>None and zmax<>None:
            zmin, zmax = min(zmin, zmax), max(zmin, zmax)
            zrange = zmin*factor, zmax*factor
        else:
            raise Exception('instance error')
        self.zrange = round(zrange[0], Layer.prec), round(zrange[1], Layer.prec)
        self.zmin = round(zrange[0], Layer.prec)
        self.zmax = round(zrange[1], Layer.prec)
        self.thickness = round(zrange[1] - zrange[0], Layer.prec)
        self.extend( [material, self.zrange] )
    def __str__(self):
        return "Layer( '%s', '%s', %s )"%(self.name, self.material, self.zrange)
    def __repr__(self):
        return "Layer( '%s', '%s', %s )"%(self.name, self.material, self.zrange)


class Layers(list):
    index = 0
    def next(self):
        layer = self[self.index]
        self.index = self.index + 1
        return layer
    def previous(self):
        self.index = self.index - 1
        layer = self[self.index]
        return layer
    def __getitem__(self, key):
        if isinstance(key, int):
            return list.__getitem__(self, key)
        for layer in self:
            if layer.name==key:
                return layer
        raise KeyError(str(key))
    def __setitem__(self, key, newlayer):
        for i, layer in enumerate(self):
            if layer.name==key or i==key:
                return list.__setitem__(self, i, newlayer)
        list.append(self, layer)
    def append(self, layer):
        if isinstance(layer, tuple):
            layer = Layer(name=layer[0], material=layer[1], zmin=min(layer[2]), zmax=max(layer[2]))
        list.append(self, layer)
       
            
            
	
		
def flatten(sequence):
    """ yield each element of an irregular list (or tuple, dict...->__instance)"""
    for el in sequence:
        if isinstance(el, (list, tuple)):
            for sub in flatten(el):
                yield sub
        else:
            yield el

class DiscretePort(object):
    index = 1
    def __init__(self, name='port', **kwargs):
        self.name = "{name}{i}".format(name=name, i=DiscretePort.index)
        DiscretePort.index += 1
        self.index = DiscretePort.index
        self.p1 = kwargs.get('P1', (0, 0, 0))
        self.p2 = kwargs.get('P2', (1, 0, 0))
    def __str__(self):
        if isinstance(self.p1, Point):
            self.p1 = self.p1.x, self.p1.y, self.p1.z
        if isinstance(self.p2, Point):
            self.p2 = self.p2.x, self.p2.y, self.p2.z
        DiscretePort.index = 1
        return "\n".join([ 'With DiscretePort',
                           '  .Reset',
                           '  .PortNumber "{i}"'.format(i=self.index-1),
                           '  .Type "SParameter"',
                           '  .Label ""',
                           '  .Impedance "50.0"',
                           '  .Voltage "1.0"',
                           '  .Current "1.0" ',
                           '  .SetP1 "False", "{p1[0]}", "{p1[1]}", "{p1[2]}" '.format(p1=self.p1),
                           '  .SetP2 "False", "{p2[0]}", "{p2[1]}", "{p2[2]}" '.format(p2=self.p2),
                           '  .InvertDirection "False" ',
                           '  .LocalCoordinates "False" ',
                           '  .Monitor "False" ',
                           '  .Radius "0.0" ',
                           '  .Create',
                           'End With',
                           ''
                        ])

class LumpedElement(object):
    index = 1
    def __init__(self, name='element', **kwargs):
        self.name = "{name}{i}".format(name=name, i=LumpedElement.index)
        LumpedElement.index += 1
        self.index = LumpedElement.index
        self.p1 = kwargs.get('P1', (0, 0, 0))
        self.p2 = kwargs.get('P2', (1, 0, 0))
        self.R = kwargs.get('R', 0)
        self.L = kwargs.get('L', 0)
        self.C = kwargs.get('C', 0)
    def __str__(self):
        return "\n".join([ 'With LumpedElement',
                           '  .Reset ',
                           '  .SetName "element1" ',
                           '  .Folder "" ',
                           '  .SetType "RLCSerial" ',
                           '  .SetR  "{r}" '.format(r=self.R),
                           '  .SetL  "{l}" '.format(l=self.L),
                           '  .SetC  "{c}" '.format(c=self.C),
                           '  .SetGs "0" ',
                           '  .SetI0 "1e-14" ',
                           '  .SetT  "300" ',
                           '  .SetP1 "False", "{p1[0]}", "{p1[1]}", "{p1[2]}" '.format(p1=self.p1),
                           '  .SetP2 "False", "{p2[0]}", "{p2[1]}", "{p2[2]}" '.format(p2=self.p2),
                           '  .SetInvert "False" ',
                           '  .SetMonitor "False" ',
                           '  .SetRadius "0.0" ',
                           '  .Create',
                           '  End With',
                           ''
                        ])


class Boundary(object):
    def __init__(self, name='boundary', **kwargs):
        self.name = "{name}".format(name=name)
        self.xmin = kwargs.get('xmin', "electric")
        self.xmax = kwargs.get('xmax', "electric")
        self.ymin = kwargs.get('ymin', "electric")
        self.ymax = kwargs.get('ymax', "electric")
        self.zmin = kwargs.get('zmin', "electric")
        self.zmax = kwargs.get('zmax', "electric")
    def __str__(self):
        return "\n".join([ 'With Boundary',
                           '  .Xmin "{self.xmin}"'.format(self=self),
                           '  .Xmax "{self.xmax}"'.format(self=self),
                           '  .Ymin "{self.ymin}"'.format(self=self),
                           '  .Ymax "{self.ymax}"'.format(self=self),
                           '  .Zmin "{self.zmin}"'.format(self=self),
                           '  .Zmax "{self.zmax}"'.format(self=self),
                           '  .Xsymmetry "none"',
                           '  .Ysymmetry "none"',
                           '  .Zsymmetry "none"',
                           '  .XminThermal "isothermal"',
                           '  .XmaxThermal "isothermal"',
                           '  .YminThermal "isothermal"',
                           '  .YmaxThermal "isothermal"',
                           '  .ZminThermal "isothermal"',
                           '  .ZmaxThermal "isothermal"',
                           '  .XsymmetryThermal "none"',
                           '  .YsymmetryThermal "none"',
                           '  .ZsymmetryThermal "none"',
                           '  .ApplyInAllDirections "False"',
                           '  .ApplyInAllDirectionsThermal "False"',
                           '  .XminTemperature ""',
                           '  .XminTemperatureType "None"',
                           '  .XmaxTemperature ""',
                           '  .XmaxTemperatureType "None"',
                           '  .YminTemperature ""',
                           '  .YminTemperatureType "None"',
                           '  .YmaxTemperature ""',
                           '  .YmaxTemperatureType "None"',
                           '  .ZminTemperature ""',
                           '  .ZminTemperatureType "None"',
                           '  .ZmaxTemperature ""',
                           '  .ZmaxTemperatureType "None"',
                           'End With',
                           ''
                           ])


class Unit(object):
    def __init__(self, name='unit', **kwargs):
        self.name = "{name}".format(name=name)
    def __str__(self):
        return "\n".join([ 'With Units',
                           '  .Geometry "um"',
                           '  .Frequency "GHz"',
                           '  .Time "ns"',
                           '  .TemperatureUnit "Kelvin"',
                           '  .Voltage "V"',
                           '  .Current "A"',
                           '  .Resistance "Ohm"',
                           '  .Conductance "Siemens"',
                           '  .Capacitance "PikoF"',
                           '  .Inductance "NanoH"',
                           '  End With',
                           '',
                           ])

            
class Material:
    tnom = 25.0
    temp = 25.0
    def __init__(self, name='material1', **kwargs):
        self.name = "{name}".format(name=name)
        self.epsilon = kwargs.get('epsilon', 1.0)
        self.kappa = kwargs.get('kappa', 0.0)
        self.type  = kwargs.get('type', "Normal")
        self.color  = kwargs.get('color', (1., 1., 0.))
        self.transparency  = kwargs.get('transparency', 0.)
        self.tc1 = kwargs.get('tc1', 0.0)
    def __str__(self):
        return "\n".join([ 'With Material',
                           '  .Reset',
                           '  .Name "{self.name}"'.format(self=self),
                           '  .FrqType "all" ',
                           '  .SetMaterialUnit "GHz", "um" ',
                           '  .Type "{self.type}" '.format(self=self),
                           '  .Epsilon "{self.epsilon}" '.format(self=self),
                           '  .Mue "1.0" ',
                           '  .Kappa "{kappa}" '.format(kappa=self.kappa/(1.0+self.tc1*(self.temp-self.tnom))),
                           '  .DispModelEps "None" ',
                           '  .DispModelMue "None" ',
                           '  .DispersiveFittingSchemeEps "General 1st" ',
                           '  .DispersiveFittingSchemeMue "General 1st" ',
                           '  .UseGeneralDispersionEps "False" ',
                           '  .UseGeneralDispersionMue "False" ',
                           '  .Colour "{self.color[0]}", "{self.color[1]}", "{self.color[2]}" '.format(self=self),
                           '  .Wireframe "False" ',
                           '  .Transparency "{self.transparency}" '.format(self=self),
                           '  .Create',
                           'End With',
			   '',
                           ])



class Mesh(object):
    def __init__(self, **kwargs):
        self.type = kwargs.get('type', 'Tetrahedral')
    def __str__(self):
        return "\n".join([ 'Mesh.MeshType "{self.type}"'.format(self=self),
                           ''
                           ])

class Solver(object):
    def __init__(self, **kwargs):
        self.frequencyrange = kwargs.get('FrequencyRange', (1, 2))
    def __str__(self):
        return "\n".join([ 'Solver.FrequencyRange "{self.frequencyrange[0]}", "{self.frequencyrange[1]}"'.format(self=self),
                           ''
                           ])


class _Extrude(list):
    _names = {'Extrude':1}
    def __init__(self, **kwargs):
        name = kwargs.get('name', 'Extrude')
        if name in Extrude._names:
            self.name = "{name}{i}".format(name=name, i=Extrude._names[name])
            Extrude._names[name] += 1
        else:
            self.name = "{name}".format(name=name)
            Extrude._names[name] = 1
        self.component = kwargs.get('component', 'component1')
        self.material = kwargs.get('material', 'PEC')
        self.zrange = kwargs.get('zrange', (0.0, 1.0))
        
    def __str__(self):
        s = []
        s.append( 'With Extrude' )
        s.append( '  .Reset' )
        s.append( '  .Name "{name}"'.format(name=self.name) )
        s.append( '  .Component "{component}"'.format(component=self.component) )
        s.append( '  .Material "{self.material}" '.format(self=self) )
        s.append( '  .Mode "Pointlist"' )
        s.append( '  .Height "{{dz{format}}}"'.format(format=CST.Formatting).format(dz=self.zrange[1]-self.zrange[0]) )
        s.append( '  .Twist "0.0"' )
        s.append( '  .Taper "0.0"' )
        s.append( '  .Origin "0.0", "0.0", "{{origin{format}}}"'.format(format=CST.Formatting).format(origin=self.zrange[0]) )
        s.append( '  .Uvector "1.0", "0.0", "0.0" ' )
        s.append( '  .Vvector "0.0", "1.0", "0.0" ' )
        if isinstance(self[0], Point):
            s.append( '  .Point "{{pt.x{format}}}", "{{pt.y{format}}}"'.format(format=CST.Formatting).format(pt=self[0]) )
        else:
            s.append( '  .Point "{{pt[0]{format}}}", "{{pt[1]{format}}}"'.format(format=CST.Formatting).format(pt=self[0]) )
        for pt in self[1:]:
            if isinstance(pt, Point):
                s.append( '  .LineTo "{{pt.x{format}}}", "{{pt.y{format}}}"'.format(format=CST.Formatting).format(pt=pt) )
            else:
                s.append( '  .LineTo "{{pt[0]{format}}}", "{{pt[1]{format}}}"'.format(format=CST.Formatting).format(pt=pt) )
        s.append( '  .Create' )
        s.append( 'End With' )
        s.append( '' )
        return "\n".join(s)

class Extrude(list):
    _names = {'Extrude':1}
    def __init__(self, **kwargs):
        name = kwargs.get('name', 'Extrude')
        if name in Extrude._names:
            self.name = "{name}{i}".format(name=name, i=Extrude._names[name])
            Extrude._names[name] += 1
        else:
            self.name = "{name}".format(name=name)
            Extrude._names[name] = 1
        self.component = kwargs.get('component', 'component1')
        self.material = kwargs.get('material', 'PEC')
        self.zrange = kwargs.get('zrange', (0.0, 1.0))
        
    def __str__(self):
        Round = lambda x: round(x*1e6)/1e6
        s = []
        s.append( 'With Extrude' )
        s.append( '  .Reset' )
        s.append( '  .Name "{name}"'.format(name=self.name) )
        s.append( '  .Component "{component}"'.format(component=self.component) )
        s.append( '  .Material "{self.material}" '.format(self=self) )
        s.append( '  .Mode "Pointlist"' )
        s.append( '  .Height "{dz}"'.format(dz=Round(self.zrange[1]-self.zrange[0])) )
        s.append( '  .Twist "0.0"' )
        s.append( '  .Taper "0.0"' )
        s.append( '  .Origin "0.0", "0.0", "{origin}"'.format(origin=Round(self.zrange[0])))
        s.append( '  .Uvector "1.0", "0.0", "0.0" ' )
        s.append( '  .Vvector "0.0", "1.0", "0.0" ' )
        pt = self[0]
        if isPoint(pt):
            x, y = Round(pt.x), Round(pt.y)
            s.append( '  .Point "{x}", "{y}"'.format(x=x, y=y) )
            x0, y0 = x, y
            for pt in self[1:]:
                x, y = Round(pt.x), Round(pt.y)
                if x==x0 and y==y0: continue
                s.append( '  .LineTo "{x}", "{y}"'.format(x=x, y=y) )
                x0, y0 = x, y
        else:
            x, y = Round(pt[0]), Round(pt[1])
            s.append( '  .Point "{x}", "{y}"'.format(x=x, y=y) )
            x0, y0 = x, y
            for pt in self[1:]:
                x, y = Round(pt[0]), Round(pt[1])
                if x==x0 and y==y0: continue
                s.append( '  .LineTo "{x}", "{y}"'.format(x=x, y=y) )
                x0, y0 = x, y
        s.append( '  .Create' )
        s.append( 'End With' )
        s.append( '' )
        return "\n".join(s)


class Group:
    def __init__(self):
        self.add = list()
        self.addItem = list()
    def Add(self, name, group):
        self.add.append((name, group))
        return self
    def AddItem(self, component, name, group):
        self.addItem.append((component, name, group))
        return self
    def __str__(self):
        s = []
        for name, group in self.add:
            s.append( 'Group.Add "{}", "{}"'.format(name, group) )
        for component, name, group in self.addItem:
            s.append( 'Group.AddItem "solid${}:{}", "{}"'.format(component, name, group) )
        s.append('')
        return "\n".join(s)

class MeshSettings:
    def __init__(self, name, size):
        self.name = name
        self.size = size
    def __str__(self):
        s = []
        s.append( 'With MeshSettings' )
        s.append( '     With .ItemMeshSettings ("group${}")'.format(self.name) )
        s.append( '          .SetMeshType "Tet"' )
        s.append( '          .Set "OctreeSizeFaces", "0"' )
        s.append( '          .Set "Size", "{}"'.format(self.size) )
        s.append( '     End With' )
        s.append( 'End With' )
        s.append('')
        return "\n".join(s)
    

try:    
    import Transform
except:
    pass
try:    
    import Solid
except:
    pass
try:    
    from Solid import Solid, Solids, Subtract, Intersect
except:
    pass
try:    
    from Generic import *
except:
    pass
