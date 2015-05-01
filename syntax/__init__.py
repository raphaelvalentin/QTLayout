from math import *
import clipper
from linalg import solve
from libarray import *
from newton import fmin as newton
from _functions import *

from CSTlib import *
from gdsii import *
import levmar
from linalg import *
from layers import *


generator = type((x for x in xrange(5)))

#__all__ = ['Point', 'Line', 'Primitive', 'Primitives', 'Vector', 'Segment', 'Path', 'Paths', 'Spline', 'linspace', 'newton']

class Point(object):
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)
        
    def __repr__(self):
        return 'Point({self.x}, {self.y})'.format(self=self)
    
    def __str__(self):
        return 'Point({self.x}, {self.y})'.format(self=self)
    
    def __add__(self, object):
        if isinstance(object, Point):
            return Point(self.x+object.x, self.y+object.y)
        if isinstance(object, (float, int)):
            return Point(self.x+object, self.y+object)
        raise Exception('%s is not a float, an integer or a Point'%repr(object))
    
    def __sub__(self, object):
        if isinstance(object, Point):
            return Point(self.x-object.x, self.y-object.y)
        if isinstance(object, (float, int)):
            return Point(self.x-object, self.y-object)
        raise Exception('%s is not a float, an integer or a Point'%repr(object))
    
    def __mul__(self, object):
        if isinstance(object, Point):
            return Point(self.x*object.x, self.y*object.y)
        if isinstance(object, (float, int)):
            return Point(self.x*object, self.y*object)
        raise Exception('%s is not a float, an integer or a Point'%repr(object))
    
    def __div__(self, object):
        if isinstance(object, Point):
            return Point(self.x/object.x, self.y/object.y)
        if isinstance(object, (float, int)):
            return Point(self.x/object, self.y/object)
        raise Exception('%s is not a float, an integer or a Point'%repr(object))
    
    def __rdiv__(self, object):
        if isinstance(object, Point):
            return Point(object.x/self.x, object.y/self.y)
        if isinstance(object, (float, int)):
            return Point(object/self.x, object/self.y)
        else:
            raise Exception('%s is not a float, an integer or a Point'%repr(object))
        
    def __pow__(self, object):
        if isinstance(object, Point):
            return Point(self.x**object.x, self.y**object.y)
        if isinstance(object, (float, int)):
            return Point(self.x**object, self.y**object)
        else:
            raise Exception('%s is not a float, an integer or a Point'%repr(object))
        
    def __radd__(self, obj):
        return self.__add__(obj)
    
    def __rsub__(self, obj):
        if isinstance(object, Point):
            return Point(object.x-self.x, object.y-self.y)
        if isinstance(object, (float, int)):
            return Point(object-self.x, object-self.y)
        else:
            raise Exception('%s is not a float, an integer or a Point'%repr(object))
        
    def __rmul__(self, obj):
        return self.__mul__(obj)
    
    def __neg__(self):
        return Point(-self.x, -self.y)
    
    def __pos__(self):
        return Point(self.x, self.y)
    
    def Rotate(self, center=None, angle=None):
        if center==None:
            center=Point(0, 0)
        if angle==None:
            angle = 0
        angle = angle%360
        if angle==0:
            return Point(self.x, self.y)
        if angle==90:
            return Point( -(self.y-center.y)+center.x,
                           (self.x-center.x)+center.y)
        if angle==180:
            return Point( -(self.x-center.x)+center.x,
                          -(self.y-center.y)+center.y)
        if angle==270:
            return Point(  (self.y-center.y)+center.x,
                          -(self.x-center.x)+center.y)
        theta = angle*pi/180.
        return Point( cos(theta)*(self.x-center.x) - sin(theta)*(self.y-center.y)+center.x,
                      sin(theta)*(self.x-center.x) + cos(theta)*(self.y-center.y)+center.y)
    
    def Mirror(self, center=None, planenormal=None):
        if center==None:
            center=Point(0, 0)
        if planenormal==None:
            planenormal=(0,0)
        if planenormal==(0,0):
            return Point(self.x, self.y)
        if planenormal==(1,0):
            return Point(center.x-self.x, self.y)
        elif planenormal==(0,1):
            return Point(self.x, center.y-self.y)
        elif planenormal==(1,1):
            return Point(center.x-self.x, center.y-self.y)
        
    def Translate(self, vector=None, angle=None, radius=None):
        if vector:
            return self + vector
        elif line and radius:
            return self + Point(radius, 0).Rotate(angle)
        return Point(self.x, self.y)
    
    def distance(self, item):
        if isinstance(item, Point):
            return sqrt((self.x-item.x)**2+(self.y-item.y)**2)
        elif item.__class__.__name__=='Line':
            vector_ortho = Point(item[0].y-item[1].y, item[1].x-item[0].x)
            point = Line(self, self+vector_ortho).Intersect(item)
            return sqrt((self.x-point.x)**2+(self.y-point.y)**2)
        
    def __eq__(self, item):
        if type(item) is Point:
            return self.x == item.x and self.y == item.y
        else:
            return False
        
    def __ne__(self, item):
        if type(item) is Point:
            return self.x <> item.x or self.y <> item.y
        else:
            return False

    def tolist(self):
        return (self.x, self.y)

    def __getitem__(self, i):
        if i==0:
            return self.x
        elif i==1:
            return self.y
        else:
            raise IndexError('list index out of range')

    def ongrid(self, grid=1e-6):
        return Point( round(self.x/grid)*grid, round(self.y/grid)*grid )


class Line(list):
    
    def __init__(self, *points):

        if len(points)==1:
            assert isinstance(points[0], generator), 'TypeError: argument is not a generator'
            points = list(points[0])
        assert len(points)==2, 'TypeError: number of points different of 2'
        assert isinstance(points[0], Point), 'TypeError: argument must be 2 Points'
        assert isinstance(points[1], Point), 'TypeError: argument must be 2 Points'
        
        assert points[0] <> points[1], 'TypeError: argument must be 2 different Points'
        assert not isinf(points[0].x) and not isinf(points[1].x), 'TypeError: arguments  must do not have infinite coordinate'
        assert not isinf(points[0].y) and not isinf(points[1].y), 'TypeError: arguments  must do not have infinite coordinate'
        
        list.__init__(self, points)

    def Rotate(self, center=None, angle=None):
            return Line(self[0].Rotate(center=center, angle=angle),
                        self[1].Rotate(center=center, angle=angle))
        
    def Translate(self, vector=None, radius=None, angle=None):
        return Line(self[0].Translate(vector=vector, radius=radius, angle=angle),
                    self[1].Translate(vector=vector, radius=radius, angle=angle))
    
    def Mirror(self, center=None, planenormal=None):
        return Line(self[0].Mirror(center=center, planenormal=planenormal),
                    self[1].Mirror(center=center, planenormal=planenormal))
    
    def __call__(self, x):
        return self.slope*x+self.intercept

    def __setitem__(self, *args, **kwargs):
        raise AttributeError("type object 'Line' does no support attribute '__setitem__'")

    def append(self, *args, **kwargs):
        raise AttributeError("type object 'Line' does no support attribute 'append'")

    def extend(self, *args, **kwargs):
        raise AttributeError("type object 'Line' does no support attribute 'extend'")

    def Intersect(self, line):
        p, q, r, s = self[0], self[1], line[0], line[1]

        det_inv = (s.x - r.x) * (q.y - p.y) - (q.x - p.x) * (s.y - r.y)
        assert det_inv<>0, 'TypeError: Lines are colinear'
        
        k1 = ((s.y - r.y) * (p.x - r.x) - (s.x - r.x) * (p.y - r.y))/det_inv
        k2 = ((q.y - p.y) * (p.x - r.x) - (q.x - p.x) * (p.y - r.y))/det_inv

        i1 = p + k1 * (q - p)
        i2 = r + k2 * (s - r)

        if (i1.x-i2.x)**2 + (i1.y-i2.y)**2 > 1e-24:
             print 'AccuracyWarning: Calculation of Intersect may be subject to inacuracy.'

        return i1

    def __str__(self):
        return '%r*x+%r'%(self.slope, self.intercept)
    
    def __repr__(self):
        return '%r*x+%r'%(self.slope, self.intercept)

    @property
    def slope(self):
        p, q = self[0], self[1]
        try:
            return (q.y - p.y) / (q.x - p.x)
        except ZeroDivisionError:
            return float('inf')
        
    @property
    def intercept(self):
        p, q = self[0], self[1]
        try:
            return (p.y * q.x - q.y * p.x) / (q.x - p.x)
        except ZeroDivisionError:
            return float('nan')

    def distance(self, obj):
        if isinstance(obj, Point):
            return obj.distance(self)
        if isinstance(obj, Line):
            assert self.slope <> obj.slope, 'TypeError: Lines are colinear. Can not calculate a distance.'
            return self[0].distance(obj)
            


class Primitive(list):
    def __init__(self, *points):
        if len(points):
            if len(points)==1:
                assert isinstance(points[0], generator), 'TypeError: argument is not a generator'
                points = list(points[0])
            if type(points[0]) == Point:
                list.__init__(self, points)
                
    def __add__(self, object):
        return Primitive(mypoint+object for mypoint in self)
    
    def __sub__(self, object):
        return Primitive(mypoint-object for mypoint in self)
    
    def __mul__(self, object):
        return Primitive(mypoint*object for mypoint in self)
    
    def __div__(self, object):
        return Primitive(mypoint/object for mypoint in self)
    
    def __rdiv__(self, object):
        return Primitive(object/mypoint for mypoint in self)
    
    def __pow__(self, object):
        return Primitive(object**mypoint for mypoint in self)
    
    def __radd__(self, obj):
        return self.__add__(obj)
    
    def __rsub__(self, obj):
        return Primitive(object-mypoint for mypoint in self)
    
    def __rmul__(self, obj):
        return self.__mul__(obj)
    
    def __neg__(self):
        return Primitive(-mypoint for mypoint in self)
    
    def __pos__(self):
        return Primitive(mypoint for mypoint in self)
    
    def Rotate(self, center=None, angle=None):
        return Primitive(mypoint.Rotate(center=center, angle=angle) for mypoint in self)
    
    def Mirror(self, center=None, planenormal=None):
        return Primitive(mypoint.Mirror(center=center, planenormal=planenormal) for mypoint in self)
    
    def Translate(self, vector=None, angle=None, radius=None):
        return Primitive(mypoint.Translate(vector=vector, angle=angle, radius=radius) for mypoint in self)
    
    def __contains__(self, item):
        if isinstance(item, Point):
            return self.PointInPolygon(self, item) in (-1, 1)
        elif isinstance(item, Primitive):
            return self.PolygonInPolygon(self, item) in (-1, 1)

    def PointInPolygon(self, point):
        return _PointInPolygon(point, self)

    def PolygonInPolygon(self, item):
        for point in self:
            if _PointInPolygon(point, item) == 0:
                return False
        return True
        return _PolygonInPolygon(self, item)
    
    @property
    def bounds(self):
        return Point(self.xmin, self.ymin), Point(self.xmax, self.ymax)
    
    def __getitem__(self, y):
        assert isinstance(y, int), 'TypeError: list indices must be integers, not str'
        return list.__getitem__(self, y%len(self))
    
    def __getslice__(self, i, j):
        return Primitive(*list.__getslice__(self, i, j))
    
    def Simplify(self, radius=1e-6):
        return Primitive(*_Simplify(self, radius))
    
    @property
    def edges(poly):
        r = []
        for i in xrange(len(poly)-1):
            r.append( Segment(poly[i], poly[i+1]) )
        r.append( Segment(poly[-1], poly[0]) )
        return r
    
    def Scale(self, center=None, scale=Point(1,1)):
        if not center:
            pt1, pt2 = self.bounds
            center = 0.5*(pt1+pt2)
        obj = Primitive()
        for point in self:
            obj.append( scale*(point-center) )
        return obj.Translate(vector=center)
    
    @property
    def xmin(self):
        return min(point.x for point in self)
    
    @property
    def xmax(self):
        return max(point.x for point in self)
    
    @property
    def ymin(self):
        return min(point.y for point in self)
    
    @property
    def ymax(self):
        return max(point.y for point in self)

    def tolist(self):
        return [p.tolist() for p in self]

    def Difference(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Difference(self, clip)
    
    def Union(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Union(self, clip)
    
    def Intersection(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Intersection(self, clip)
    
    def Xor(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Xor(self, clip)

    def __eq__(self, item):
        if type(item) is Primitive:
            for point1, point2 in zip(self, item):
                if point1<>point2:
                    return False
            return True
        else:
            return False
            
    def __ne__(self, item):
        if type(item) is Primitive:
            for point1, point2 in zip(self, item):
                if point1==point2:
                    return False
        else:
            return False
        
    def ongrid(self, grid):
        return Primitive(mypoint.ongrid(grid) for mypoint in self)

class Primitives(list):
    def __init__(self, *primitives):
        if len(primitives):
            if len(primitives)==1:
                assert isinstance(primitives[0], generator), 'TypeError: argument is not a generator'
                primitives = list(primitives[0])
            if type(primitives[0]) == Primitive:
                list.__init__(self, primitives)

        for primitive in primitives:
            if not(isinstance(primitive, Primitive)):
                raise Exception('at least, one element is not a primitive')
        list.__init__(self, list(primitives))
        
    def Mirror(self, center=None, planenormal=None):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive.Mirror(center=center, planenormal=planenormal) )
        return obj
    
    def Translate(self, vector=None, radius=None, angle=None):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive.Translate(vector=vector, radius=radius, angle=angle) )
        return obj
    
    def Rotate(self, center=None, angle=None):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive.Rotate(center=center, angle=angle) )
        return obj

    def Scale(self, center=None, scale=Point(1,1)):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive.Scale(center=center, scale=scale) )
        return obj
    
    @property
    def bounds(self):
        return Point(self.xmin, self.ymin), Point(self.xmax, self.ymax)
    
    def __mul__(self, object):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive*object )
        return obj
    
    def __rmul__(self, object):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive*object )
        return obj
    
    def __div__(self, object):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive/object )
        return obj
    
    def __rdiv__(self, object):
        obj = Primitives()
        for primitive in self:
            obj.append( object/primitive )
        return obj
    
    def __add__(self, object):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive+object )
        return obj

    def __radd__(self, obj):
        return self.__add__(obj)
    
    def __sub__(self, object):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive-object )
        return obj
    
    @property
    def xmin(self):
        return min(primitive.xmin for primitive in self)
    
    @property
    def xmax(self):
        return max(primitive.xmax for primitive in self)
    
    @property
    def ymin(self):
        return min(primitive.ymin for primitive in self)
    
    @property
    def ymax(self):
        return max(primitive.ymax for primitive in self)

    def Simplify(self, radius=1e-6):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive.Simplify(radius=radius) )
        return obj

    def Difference(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Difference(self, clip)
    
    def Union(self, *objs):
        clip = Primitives()
        if len(objs)>0:
            for obj in objs:
                if isinstance(obj, Primitive):
                    clip.append(obj)
                elif isinstance(obj, Primitives):
                    for primitive in obj:
                        clip.append(primitive)
            return Union(self, clip)
        else:
            return Union(self[0], Primitives(*self[1:]))
    
    def Intersection(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Intersection(self, clip)
    
    def Xor(self, *objs):
        clip = Primitives()
        for obj in objs:
            if isinstance(obj, Primitive):
                clip.append(obj)
            elif isinstance(obj, Primitives):
                for primitive in obj:
                    clip.append(primitive)
        return Xor(self, clip)

    def ongrid(self, grid):
        obj = Primitives()
        for primitive in self:
            obj.append( primitive.ongrid(grid) )
        return obj
        
        
# clipper wrap
def Intersection(subject, clip):
    return _Clip(subject, clip, clipper.ClipType.Intersection)
def Union(subject, clip):
    return _Clip(subject, clip, clipper.ClipType.Union)
def Difference(subject, clip):
    return _Clip(subject, clip, clipper.ClipType.Difference)
def Xor(subject, clip):
    return _Clip(subject, clip, clipper.ClipType.Xor)
    
def _Clip(subject, clip, cliptype):
    c = clipper.Clipper()
    if isinstance(subject, Primitive):
        polygon = [clipper.Point(point.x, point.y) for point in subject]
        c.AddPolygon(polygon, clipper.PolyType.Subject)
    elif isinstance(subject, Primitives):
        for primitive in subject:
            polygon = [clipper.Point(point.x, point.y) for point in primitive]
            c.AddPolygon(polygon, clipper.PolyType.Subject)
    if isinstance(clip, Primitive):
        polygon = [clipper.Point(point.x, point.y) for point in clip]
        c.AddPolygon(polygon, clipper.PolyType.Subject)
    elif isinstance(clip, Primitives):
        for primitive in clip:
            polygon = [clipper.Point(point.x, point.y) for point in primitive]
            c.AddPolygon(polygon, clipper.PolyType.Clip)
    solution = []
    pft = clipper.PolyFillType.NonZero
    result = c.Execute(cliptype, solution, pft, pft)
    primitives = Primitives()
    for bloc in solution:
        p = Primitive()
        for point in bloc:
            p.append( Point(float(point.x), float(point.y)) )
        primitives.append(p)
    return primitives
                        

class Vector(Point):
    def orthogonal(self):   
        # +90deg
        return Vector(-self.y, self.x)

    def norm(self):
        return sqrt(self.x*self.x + self.y*self.y)

    def __repr__(self):
        return 'Vector({self.x}, {self.y})'.format(self=self)
    def __str__(self):
        return 'Vector({self.x}, {self.y})'.format(self=self)


class Segment(list):
    def __init__(self, *points):
        assert len(points) == 2, 'TypeError: Segment() takes exactly 2 arguments'
        #for point in points:
        #    if not(isinstance(point, Point)):
        #        raise Exception('at least, one element is not a Point')
        list.__init__(self, list(points))
    def Intersect(self, segment):
        """  The main function that returns True if line segment 'p1q1'
             and 'p2q2' intersect.
        """
        
        def onSegment(p, q, r):
            """ Given three colinear pos p, q, r, the function checks if
            point q lies on line segment 'pr'
            """
            if (q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)):
                return True
            return False

        def orientation(p, q, r):
            """ To find orientation of ordered triplet (p, q, r).
             The function returns following values
             0 --> p, q and r are colinear
             1 --> Clockwise
             2 --> Counterclockwise
             """
            #See 10th slides from following link for derivation of the formula
            #http:#www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
            val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
            if val == 0:
                return 0  # colinear
            if val > 0: # clock or counterclock wise
                return 1
            else:
                return 2
            
        p1, q1, p2, q2 = self[0], self[1], segment[0], segment[1]
        # Find the four orientations needed for general and
        # special cases
        o1 = orientation(p1, q1, p2)
        o2 = orientation(p1, q1, q2)
        o3 = orientation(p2, q2, p1)
        o4 = orientation(p2, q2, q1)
        # General case
        if (o1 != o2 and o3 != o4):
            return Line(p1, q1).Intersect( Line(p2, q2) )
        # Special Cases
        # p1, q1 and p2 are colinear and p2 lies on segment p1q1
        if (o1 == 0 and onSegment(p1, p2, q1)):
            return p2
        # p1, q1 and p2 are colinear and q2 lies on segment p1q1
        if (o2 == 0 and onSegment(p1, q2, q1)):
            return q2
        # p2, q2 and p1 are colinear and p1 lies on segment p2q2
        if (o3 == 0 and onSegment(p2, p1, q2)):
            return p1
         # p2, q2 and q1 are colinear and q1 lies on segment p2q2
        if (o4 == 0 and onSegment(p2, q1, q2)):
            return q1
        return None # Doesn't fall in any of the above cases
        
    def Rotate(self, center=None, angle=None):
            return Segment(self[0].Rotate(center=center, angle=angle),
                           self[1].Rotate(center=center, angle=angle))
        
    def Translate(self, vector=None, radius=None, angle=None):
        return Segment(self[0].Translate(vector=vector, radius=radius, angle=angle),
                       self[1].Translate(vector=vector, radius=radius, angle=angle))
    
    def Mirror(self, center=None, planenormal=None):
        return Segment(self[0].Mirror(center=center, planenormal=planenormal),
                       self[1].Mirror(center=center, planenormal=planenormal))
    
    def Length(self):
        return self[0].distance(self[1])


class Path(list):
    def __init__(self, *points):
        list.__init__(self, list(points))
        
    @property
    def edges(poly):
        r = []
        for i in xrange(len(poly)-1):
            r.append( Segment(poly[i], poly[i+1]) )
        return r
    
    @property
    def vectors(poly):
        r = []
        for i in xrange(len(poly)-1):
            r.append( Vector(poly[i+1].x-poly[i].x, poly[i+1].y-poly[i].y) )
        return r
    
    def Simplify(self, radius=1e-6):
        return Path(*_Simplify(self, radius))
    
    def enlarge(self, width):
        primitive = Primitive()

        vector = self.vectors[0].orthogonal()
        primitive.append( self[0]-0.5*width*vector/vector.norm() )

        vector = self.vectors[0].orthogonal()
        primitive.append( self[0]+0.5*width*vector/vector.norm() )
    
        for i in xrange(len(self)-2):
            vector = self.vectors[i].orthogonal()
            line3 = Line(self[i], self[i+1]).Translate(vector=0.5*width*vector/vector.norm())
            vector = self.vectors[i+1].orthogonal()
            line4 = Line(self[i+1], self[i+2]).Translate(vector=0.5*width*vector/vector.norm())
            try:
                primitive.append( line3.Intersect(line4) )
            except:
                # colinear
                pass

        vector = self.vectors[-1].orthogonal()
        primitive.append( self[-1]+0.5*width*vector/vector.norm() )
    
        vector = self.vectors[-1].orthogonal()
        primitive.append( self[-1]-0.5*width*vector/vector.norm() )
   
        for i in xrange(len(self)-3, -1, -1):
            vector = self.vectors[i].orthogonal()
            line3 = Line(self[i], self[i+1]).Translate(vector=-0.5*width*vector/vector.norm())
            vector = self.vectors[i+1].orthogonal()
            line4 = Line(self[i+1], self[i+2]).Translate(vector=-0.5*width*vector/vector.norm())
            try:
                primitive.append( line3.Intersect(line4) )
            except:
                # colinear
                pass

        return primitive
    
    def Rotate(self, center=None, angle=None):
        return Path(*[mypoint.Rotate(center=center, angle=angle) for mypoint in self])
    
    def Mirror(self, center=None, planenormal=None):
        return Path(*[mypoint.Mirror(center=center, planenormal=planenormal) for mypoint in self])
    
    def Translate(self, vector=None):
        return Path(*[mypoint.Translate(vector=vector) for mypoint in self])
    
    def __add__(self, object):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[mypoint+object for mypoint in self])
    
    def __sub__(self, object):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[mypoint-object for mypoint in self])
    
    def __mul__(self, object):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[mypoint*object for mypoint in self])
    
    def __div__(self, object):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[mypoint/object for mypoint in self])
    
    def __rdiv__(self, object):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[object/mypoint for mypoint in self])
    
    def __pow__(self, object):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[object**mypoint for mypoint in self])
    
    def __radd__(self, obj):
        return self.__add__(obj)
    
    def __rsub__(self, obj):
        assert isinstance(object, Point) , 'TypeError: %s is not a Point'%repr(object)
        return Path(*[object-mypoint for mypoint in self])
    
    def __rmul__(self, obj):
        return self.__mul__(obj)
    
    def __neg__(self):
        return Path(*[-mypoint for mypoint in self])
    
    def __pos__(self):
        return Path(*[mypoint for mypoint in self])
    
    @property
    def length(self):
        d = 0.0
        point0 = self[0]
        for point1 in self[1:]:
            d += point0.distance(point1)
            point0 = point1
        return d
    
    @property
    def xmin(self):
        return min(point.x for point in self)
    
    @property
    def xmax(self):
        return max(point.x for point in self)
    
    @property
    def ymin(self):
        return min(point.y for point in self)
    
    @property
    def ymax(self):
        return max(point.y for point in self)
    
class Paths(list):
    def __init__(self, *paths):
        for path in paths:
            if not(isinstance(path, Path)):
                raise Exception('at least, one element is not a Path')
        list.__init__(self, list(paths))
        
def polyval(p,x):
    y=0.
    l=len(p)
    for i,v in enumerate(p):
        y+=v*x**(l-i-1)
    return y

def Spline(curve1, curve2, n=100):
    n = int(n)
    if isinstance(curve1, Segment) and isinstance(curve2, Segment): 
        x0, y0 = curve1[1].x, curve1[1].y
        x1, y1 = curve2[0].x, curve2[0].y
        slope0 = (curve1[1].y-curve1[0].y)/(curve1[1].x-curve1[0].x)
        slope1 = (curve2[1].y-curve2[0].y)/(curve2[1].x-curve2[0].x)
        A = array([[x0**3, x0**2, x0, 1.0],
                   [x1**3, x1**2, x1, 1.0],
                   [3.0*x0**2, 2.0*x0, 1.0, 0.0],
                   [3.0*x1**2, 2.0*x1, 1.0, 0.0]])
        B = array([y0, y1, slope0, slope1])
        X = solve(A, B)
        curve3 = []
        for x in linspace(x0, x1, n):
            curve3.append( Point(x, polyval(X, x)) )
        return Path(*curve3)


####### functions

def linspace(start, stop, num=0, step=None):
    start = float(start)
    stop = float(stop)
    if step == None:
        if num == 1:
            step = stop-start
        else:
            step = (stop-start)/(float(num)-1.0)
    else:
        step = float(step)
        if step == 0.0:
            num = 0
        else:
            num = (stop-start)/step+1.0
    return array(trunc([i*step+start for i in xrange(int(num + pow(2,-45)))]))

def trunc(x, ndigits=15):
    """ trunc(number/sequence, ndigits=15) -> float/list point number(s) """
    def __trunc__(x):
        if x<>0:
            return round(x, -int(round(log10(abs(x)))+1)+ndigits)
	else:
	    return 0.0
    if hasattr(x, "__iter__"):
        return 	map(__trunc__, x)
    else: 
        return __trunc__(x)

