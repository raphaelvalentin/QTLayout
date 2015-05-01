######################################
# fonctions primitive
fmin = newton.NewtonRaphson

execfile('./Projects/primitives.py') 

def Bridge1x2(width1=12.0, width2=10.0, spacing=3.0):
    
    paths = []

    def f(dx, paths=None):

        y0l = width1+0.5*width2+spacing
        y1l = y0l-width2
        y2l = y1l-spacing
        y3l = y2l-width1
        y4l = y3l-spacing
        y5l = y4l-width1

        y0r = width1+0.5*width2+spacing
        y1r = y0r-width1
        y2r = y1r-spacing
        y3r = y2r-width1
        y4r = y3r-spacing
        y5r = y4r-width2

        path1 = Path()
        line1 = Line( Point(-100, 0.5*(y0l+y1l), 0), Point(100, 0.5*(y0l+y1l), 0) )
        line2 = Line( Point(-100, 100, 0), Point(100, -100, 0) )
        point1 = line1.Intercept(line2)
        line1 = Line( Point(-100, 0.5*(y4r+y5r), 0), Point(100, 0.5*(y4r+y5r), 0) )
        point2 = line1.Intercept(line2)

        path1.append( Point(-100, point1.y, 0) )
        path1.append(point1)
        path1.append(point2)
        path1.append( Point(100, point2.y, 0) )
        primitive1 = path1.enlarge(width=width2)
        if isinstance(paths, list):
            paths.append( path1 )

        path2 = Path()
        point3 = Point(point1.x, 0.5*(y2l+y3l), 0)
        point4 = Point(point2.x, 0.5*(y0r+y1r), 0)
        path2.append(Point(-100, point3.y, 0))
        path2.append(point3-Point(dx, 0, 0))
        path2.append(point4-Point(dx, 0, 0))
        path2.append(Point(100, point4.y, 0))
        primitive2 = path2.enlarge(width=width1)
        if isinstance(paths, list):
            paths.append( path2 )

        path3 = Path()
        point5 = Point(point1.x, 0.5*(y4l+y5l), 0)
        point6 = Point(point2.x, 0.5*(y2r+y3r), 0)
        path3.append(Point(-100, point5.y, 0))
        path3.append(point5+Point(dx, 0, 0))
        path3.append(point6+Point(dx, 0, 0))
        path3.append(Point(100, point6.y, 0))
        primitive3 = path3.enlarge(width=width1)
        if isinstance(paths, list):
            paths.append( path3 )

        line1 = Line( Point(-100, 100, 0), Point(100, -100, 0) )
        line2 = Line( primitive2[6], primitive2[7] )
        line3 = Line( primitive3[2], primitive3[3] )

        point1 = line1.Intercept(line2)
        point2 = line1.Intercept(line3)
        return point1.distance(point2)-spacing

    dx = fmin(f=f, x0=spacing, verbose=False).run()
    f(dx, paths)

    # adjust port points
    xmin = min(paths[0][1].x, paths[1][1].x, paths[2][1].x)
    paths[0][0].x = xmin-3.5
    paths[1][0].x = xmin-3.5
    paths[2][0].x = xmin-3.5

    xmax = max(paths[0][-2].x, paths[1][-2].x, paths[2][-2].x)
    paths[0][-1].x = xmax+3.5
    paths[1][-1].x = xmax+3.5
    paths[2][-1].x = xmax+3.5
    
    return paths[0], paths[1], paths[2]


def Bridge1x1(width=12.0, spacing=3.0):
    line1 = Line( Point(-100, 0.5*(spacing+width), 0), Point(100, 0.5*(spacing+width), 0) )
    line2 = Line( Point(-100, 100, 0), Point(100, -100, 0) )
    point1 = line1.Intercept(line2)
    point2 = point1.Mirror(planenormal=(1,1,0))

    path1 = Path()
    path1.append( point1-Point(100, 1e-12, 0) )
    path1.append( point1 )
    path1.append( point2 )
    path1.append( point2+Point(100, 1e-12, 0) )

    obj1 = path1.enlarge(width=width)
    x1, x2 = obj1[2].x, obj1[-1].x
    dx = 0.5*abs(x1-x2)
    
    path1 = Path()
    path1.append( point1-Point(dx+3.5, 1e-12, 0) )
    path1.append( point1 )
    path1.append( point2 )
    path1.append( point2+Point(dx+3.5, 1e-12, 0) )
    
    return path1, path1.Mirror(planenormal=(0,1,0))

        


 



"""
path1, path2, path3 = Bridge1x2()

a = path1.enlarge(width=width2)
b = path2.enlarge(width=width1)
c = path3.enlarge(width=width1)
layer1 = Primitives()
layer1.append(a)
layer1.append(b)
layer1.append(c)
"""

