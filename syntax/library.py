from primitives import *

#------------------------------------------------------------------------------
# Miscellaneous global functions
#------------------------------------------------------------------------------

def SlopesEqual(pt1, pt2, pt3, pt4):
    return (pt1.y-pt2.y)*(pt3.x-pt4.x) == (pt1.x-pt2.x)*(pt3.y-pt4.y)

def IsHorizontal(pt1, pt2):
    return pt1.y==pt2.y

def IsVertical(pt1, pt2):
    return pt1.x==pt2.x


def PointInPolygon_old(self, item):
    cn = 0    # the  crossing number counter
    # loop through all edges of the polygon
    for i in xrange(len(self)-1):  # edge from V[i]  to V[i+1]
        if (((self[i].y <= item.y) and (self[i+1].y > item.y)) or ((self[i].y > item.y) and (self[i+1].y <= item.y))):
            # compute  the actual edge-ray intersect x-coordinate
            vt = (item.y  - self[i].y) / (self[i+1].y - self[i].y)
            if (item.x <  self[i].x + vt * (self[i+1].x - self[i].x)) : # P.x < intersect
                cn = cn + 1   # a valid crossing of y=P.y right of P.x
    return bool(cn%2)  # 0 if even (out), and 1 if  odd (in)



def PointInPolygon(pt, path):
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


def PolygonInPolygon(polygon1, polygon2):
    for point in polygon1:
        result = PointInPolygon(point, polygon2)
        if result in (-1, 1):
            return result
    return 0            


# clip
from clip import *

def difference(primitive1, primitive2):
    a = []
    for p in primitive1:
        a.append( (p.x, p.y) )
    b = []
    for p in primitive2:
        b.append( (p.x, p.y) )
    c = clip_polygon(a, b, 'difference')
    d = Primitive()
    for p in c:
        for point in p.points:
             d.append( Point(point[0], point[1], 0) ) 
    return d

def union(primitive1, primitive2):
    a = []
    for p in primitive1:
        a.append( (p.x, p.y) )
    b = []
    for p in primitive2:
        b.append( (p.x, p.y) )
    c = clip_polygon(a, b, 'union')
    d = Primitive()
    for p in c:
        for point in p.points:
             d.append( Point(point[0], point[1], 0) ) 
    return d
    

def intersection(primitive1, primitive2):
    a = []
    for p in primitive1:
        a.append( (p.x, p.y) )
    b = []
    for p in primitive2:
        b.append( (p.x, p.y) )
    c = clip_polygon(a, b, 'intersection')
    d = Primitive()
    for p in c:
        for point in p.points:
             d.append( Point(point[0], point[1], 0) ) 
    return d

def reversed_diff(primitive1, primitive2):
    a = []
    for p in primitive1:
        a.append( (p.x, p.y) )
    b = []
    for p in primitive2:
        b.append( (p.x, p.y) )
    c = clip_polygon(a, b, 'reversed-diff')
    d = Primitive()
    for p in c:
        for point in p.points:
             d.append( Point(point[0], point[1], 0) ) 
    return d



