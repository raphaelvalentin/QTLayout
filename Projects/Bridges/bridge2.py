

def Bridge1X1( width1 = 12.0, width2 = 6.0, spacing = 2.0, angle = 45.0):

    # trace1
    line1 = Line( Point(-100, width1+0.5*spacing, 0), Point(100, width1+0.5*spacing, 0) )
    line2 = Line( Point(-100, 0.5*spacing, 0), Point(100, 0.5*spacing, 0) )

    # trace2
    line3 = Line( Point(-100, 0.25*(width1+width2), 0), Point(100, 0.25*(width1+width2), 0) )
    line4 = Line( Point(-100, -0.25*(width1+width2), 0), Point(100, -0.25*(width1+width2), 0) )
    line3 = line3.Rotate(angle=(0, 0, -angle))
    line4 = line4.Rotate(angle=(0, 0, -angle))

    # trace3
    line5 = Line( Point(-100, -0.5*spacing, 0), Point(100, -0.5*spacing, 0) )
    line6 = Line( Point(-100, -width2-0.5*spacing, 0), Point(100, -width2-0.5*spacing, 0) )

    #intersects
    point1 = line1.Intercept(line3)
    point2 = line3.Intercept(line5)
    point3 = line4.Intercept(line6)
    point4 = line2.Intercept(line4)

    # create bridge
    bridge = Primitive()
    bridge.append( point4 )
    bridge.append( Point(point4.x, point1.y, 0) )
    bridge.append( point1 )
    bridge.append( point2 )
    bridge.append( Point(point2.x, point3.y, 0) )
    bridge.append( point3 )
    return bridge




#layer1 = Primitives()
#layer1.append(Bridge1X1(width1=12, width2=12))



