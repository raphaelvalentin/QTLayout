
def Bridge2x2(width=(10, 8, 9, 7), spacing=2, angle=45):

    width1, width2, width3, width4 = width

    line1 = Line(Point(-100, width1+width2+1.5*spacing), Point(100,width1+width2+1.5*spacing))
    line2 = Line(Point(-100, width2+1.5*spacing), Point(100, width2+1.5*spacing))
    line3 = Line(Point(-100, 0.5*(width1+width3)+0.5*spacing), Point(100,0.5*(width1+width3)+0.5*spacing)).Rotate(angle=-angle)
    line4 = Line(Point(-100, 0.5*spacing), Point(100, 0.5*spacing)).Rotate(angle=-angle)
    line5 = Line(Point(-100, -0.5*spacing), Point(100, -0.5*spacing))
    line6 = Line(Point(-100, -0.5*spacing-width3), Point(100, -0.5*spacing-width3))

    point1 = line1.Intersect(line3)
    point2 = line3.Intersect(line5)
    point3 = line4.Intersect(line6)
    point4 = line2.Intersect(line4)

    bridge1 = Primitive( )
    bridge1.append( point4 )
    bridge1.append( Point(point4.x-3, point4.y) )
    bridge1.append( Point(point4.x-3, point1.y) )
    bridge1.append( point1 )
    bridge1.append( point2 )
    bridge1.append( Point(point2.x+3, point2.y) )
    bridge1.append( Point(point2.x+3, point3.y) )
    bridge1.append( point3 )    

    line1 = Line(Point(-100, width2+0.5*spacing), Point(100, width2+0.5*spacing))
    line2 = Line(Point(-100, 0.5*spacing), Point(100, 0.5*spacing))
    line3 = Line(Point(-100, -0.5*spacing), Point(100, -0.5*spacing)).Rotate(angle=-angle)
    line4 = Line(Point(-100, -0.5*(width2+width4)-0.5*spacing), Point(100, -0.5*(width2+width4)-0.5*spacing)).Rotate(angle=-angle)
    line5 = Line(Point(-100, -width3-1.5*spacing), Point(100, -width3-1.5*spacing))
    line6 = Line(Point(-100, -width3-width4-1.5*spacing), Point(100, -width3-width4-1.5*spacing))

    point1 = line1.Intersect(line3)
    point2 = line3.Intersect(line5)
    point3 = line4.Intersect(line6)
    point4 = line2.Intersect(line4)

    bridge2 = Primitive( )
    bridge2.append( point4 )
    bridge2.append( Point(point4.x-3, point4.y) )
    bridge2.append( Point(point4.x-3, point1.y) )
    bridge2.append( point1 )
    bridge2.append( point2 )
    bridge2.append( Point(point2.x+3, point2.y) )
    bridge2.append( Point(point2.x+3, point3.y) )
    bridge2.append( point3 )

    x1 = max(pt.x for pt in bridge1)
    x2 = min(pt.x for pt in bridge2)

    bridge1 = bridge1.Translate(vector=-0.5*Point(x1+x2, 0))
    bridge2 = bridge2.Translate(vector=-0.5*Point(x1+x2, 0))
    return bridge1, bridge2


#layer1 = Primitives()
#layer2 = Primitives()
#bridge1, bridge2 = Bridge2x2()
#layer1.append(bridge1)
#layer1.append(bridge2)
#layer2.append( bridge1.Mirror(planenormal=(1,0,0)) )
#layer2.append( bridge2.Mirror(planenormal=(1,0,0)) )



def Bridge1x1(width=(10, 8), spacing=2, angle=45):

    if isinstance(width, (float, int)):
        width1, width2, width3 = width, width, width
    elif isinstance(width, (list, tuple)):
        if len(width)==1:
            width1, width2, width3 = width[0], width[0], width[0]
        elif len(width)==2:
            width1, width2 = width
            width3 = 0.5*(width1+width2)
        elif len(width)==3:
            width1, width3, width2  = width


    line1 = Line(Point(-100, width1+0.5*spacing), Point(100,width1+0.5*spacing))
    line2 = Line(Point(-100, 0.5*spacing), Point(100, 0.5*spacing))
    line3 = Line(Point(-100, 0.5*width3), Point(100,0.5*width3)).Rotate(angle=-angle)
    line4 = Line(Point(-100, -0.5*width3), Point(100, -0.5*width3)).Rotate(angle=-angle)
    line5 = Line(Point(-100, -0.5*spacing), Point(100, -0.5*spacing))
    line6 = Line(Point(-100, -0.5*spacing-width2), Point(100, -0.5*spacing-width2))

    point1 = line1.Intersect(line3)
    point2 = line3.Intersect(line5)
    point3 = line4.Intersect(line6)
    point4 = line2.Intersect(line4)

    bridge1 = Primitive( )
    bridge1.append( point4 )
    bridge1.append( Point(point4.x-3, point4.y) )
    bridge1.append( Point(point4.x-3, point1.y) )
    bridge1.append( point1 )
    bridge1.append( point2 )
    bridge1.append( Point(point2.x+3, point2.y) )
    bridge1.append( Point(point2.x+3, point3.y) )
    bridge1.append( point3 )    
    return bridge1


#layer1 = Primitives()
#bridge1, bridge2 = Bridge2x2()
#layer1.append(bridge1)
#layer1.append(bridge2)

#layer2 = Primitives()
#layer2.append(bridge1.Mirror(planenormal=(1,0,0)))
#layer2.append(bridge2.Mirror(planenormal=(1,0,0)))



