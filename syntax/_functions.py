from math import *

__all__ = ['_PointInPolygon', '_PolygonInPolygon', '_Simplify', '_Clockwise', '_Distance']

def _PointInPolygon(pt, path):
    #returns 0 if false, +1 if true, -1 if pt ON polygon boundary
    #http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.88.5498&rep=rep1&type=pdf
    result = 0
    cnt = len(path)
    if cnt < 3:
        return 0
    ip = path[0]
    for i in xrange(1, cnt+1, 1):
        ipNext = path[0] if i == cnt else path[i]
        if ipNext.y == pt.y:
            if (ipNext.x == pt.x) or (ip.y == pt.y and ((ipNext.x > pt.x) == (ip.x < pt.x))):
                return -1
        if (ip.y < pt.y) != (ipNext.y < pt.y):
            if ip.x >= pt.x:
                if ipNext.x > pt.x:
                    result = 1 - result
                else:
                    d = (ip.x - pt.x) * (ipNext.y - pt.y) - (ipNext.x - pt.x) * (ip.y - pt.y)
                    if not(d):
                        return -1;
                    if (d > 0) == (ipNext.y > ip.y):
                        result = 1 - result
            else:
                if (ipNext.x > pt.x):
                    d = (ip.x - pt.x) * (ipNext.y - pt.y) - (ipNext.x - pt.x) * (ip.y - pt.y)
                    if not(d):
                        return -1
                    if (d > 0) == (ipNext.y > ip.y):
                        result = 1 - result
        ip = ipNext
    return result

def _PointInPolygon_old(self, item):
    cn = 0    # the  crossing number counter
    # loop through all edges of the polygon
    for i in xrange(len(self)-1):  # edge from V[i]  to V[i+1]
        if (((self[i].y <= item.y) and (self[i+1].y > item.y)) or ((self[i].y > item.y) and (self[i+1].y <= item.y))):
            # compute  the actual edge-ray intersect x-coordinate
            vt = (item.y  - self[i].y) / (self[i+1].y - self[i].y)
            if (item.x <  self[i].x + vt * (self[i+1].x - self[i].x)) : # P.x < intersect
                cn = cn + 1   # a valid crossing of y=P.y right of P.x
    return bool(cn%2)  # 0 if even (out), and 1 if  odd (in)

def _Intersect(s1, s2, c1, c2):
    """Test the intersection between two lines (two pairs of coordinates for two points).
    Return the coordinates for the intersection and the subject and clipper alphas if the test passes.
    Algorithm based on: http://paulbourke.net/geometry/lineline2d/
    """

    den = (c2.y - c1.y) * (s2.x - s1.x) - (c2.x - c1.x) * (s2.y - s1.y)
    if den == 0:
        raise TypeError('TypeError: Lines are colinear')

    us = ((c2.x - c1.x) * (s1.y - c1.y) - (c2.y - c1.y) * (s1.x - c1.x)) / den
    uc = ((s2.x - s1.x) * (s1.y - c1.y) - (s2.y - s1.y) * (s1.x - c1.x)) / den

    if (0 <= us <= 1) and (0 <= uc <= 1):
        x = s1.x + us * (s2.x - s1.x)
        y = s1.y + us * (s2.y - s1.y)
        return (x, y)
    raise TypeError('TypeError: No Intersection found')
    
def _PolygonInPolygon(polygon1, polygon2):
    for point in polygon1:
        if _PointInPolygon(point, polygon2) == 0:
            return False
    return True

def _PolygonOutPolygon(polygon1, polygon2):
    for point in polygon1:
        if _PointInPolygon(point, polygon2) in (-1, 1):
            return False
    for point in polygon2:
        if _PointInPolygon(point, polygon1) in (-1, 1):
            return False
    segments1 = []
    for i in xrange(len(polygon1)-1):
        segments1.append( (polygon1[i], polygon1[i+1]) )
    segments1.append( (polygon1[-1], polygon1[0]) )
    segments2 = []
    for i in xrange(len(polygon2)-1):
        segments2.append( (polygon2[i], polygon2[i+1]) )
    segments2.append( (polygon2[-1], polygon2[0]) )

    for segment1 in segments1:
        for segment2 in segments2:
            try:
                _Intersect(segment1[0], segment1[1], segment2[0], segment2[1])
                return False
            except:
                pass
    return True

def _Simplify(points, radius=1e-6):
    shape = []
    i=0
    while i<len(points)-1:
        if points[i].distance(points[i+1])>radius:
            shape.append(points[i])
            i = i+1
        else:
            shape.append(0.5*(points[i]+points[i+1]))
            i = i+2
    if self[-1].distance(shape[-1])>radius:
        shape.append(points[-1])
    return shape

def _Distance(point1, point2):
    return sqrt( (point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 )

def _Clockwise(point1, point2, point3):
    """ A simple test of vertex ordering for convex polygons
        is based on considerations of the cross product between adjacent edges.
        If the cross product is positive then it rises above the plane
        (z axis up out of the plane) and if negative then the cross product is into the plane. """
    z = (point2[0]-point1[0])*(point3[1]-point2[1])-(point2[1]-point1[1])*(point3[0]-point2[0])
    return z<0

