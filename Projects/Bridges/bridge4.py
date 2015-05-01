
def Bridge3x3(width = (6.5, 6.5, 11.0, 6.5, 6.5, 11.0), spacing = 2.01/1.1, angle = 45.0):
    width1, width2, width3, width4, width5, width6 = width
    #
    y0 = 0.0
    y1 = y0 + 0.5*spacing
    y2 = y1 + width3
    y3 = y2 + spacing
    y4 = y3 + width2
    y5 = y4 + spacing
    y6 = y5 + width1
    #
    y_1 = y0 - 0.5*spacing
    y_2 = y_1 - width4
    y_3 = y_2 - spacing
    y_4 = y_3 - width5
    y_5 = y_4 - spacing
    y_6 = y_5 - width6
    #
    width7 = 0.5*(width1+width4)
    width8 = 0.5*(width2+width5)
    width9 = 0.5*(width6+width6)
    #
    tw1 = width7+width8+width9+2*spacing
    #
    line1 = Line(Point(-100, y6), Point(100, y6))
    line2 = Line(Point(-100, y5), Point(100, y5))
    line3 = Line(Point(-100, 0.5*tw1), Point(100,0.5*tw1)).Rotate(angle=-angle)
    line4 = Line(Point(-100, 0.5*tw1-width7), Point(100, 0.5*tw1-width7)).Rotate(angle=-angle)
    line5 = Line(Point(-100, y_1), Point(100, y_1))
    line6 = Line(Point(-100, y_2), Point(100, y_2))
    #
    point1 = line1.Intersect(line3)
    point2 = line3.Intersect(line5)
    point3 = line4.Intersect(line6)
    point4 = line2.Intersect(line4)
    #
    bridge1 = Primitive( )
    bridge1.append( point4 )
    bridge1.append( Point(point4.x-3, point4.y) )
    bridge1.append( Point(point4.x-3, point1.y) )
    bridge1.append( point1 )
    bridge1.append( point2 )
    bridge1.append( Point(point2.x+3, point2.y) )
    bridge1.append( Point(point2.x+3, point3.y) )
    bridge1.append( point3 )    
    
    #
    line1 = Line(Point(-100, y4), Point(100, y4))
    line2 = Line(Point(-100, y3), Point(100, y3))
    line3 = Line(Point(-100, 0.5*tw1-width7-spacing), Point(100, 0.5*tw1-width7-spacing)).Rotate(angle=-angle)
    line4 = Line(Point(-100, 0.5*tw1-width7-width8-spacing), Point(100, 0.5*tw1-width7-width8-spacing)).Rotate(angle=-angle)
    line5 = Line(Point(-100, y_3), Point(100, y_3))
    line6 = Line(Point(-100, y_4), Point(100, y_4))
    #
    point1 = line1.Intersect(line3)
    point2 = line3.Intersect(line5)
    point3 = line4.Intersect(line6)
    point4 = line2.Intersect(line4)
    #
    bridge2 = Primitive( )
    bridge2.append( point4 )
    bridge2.append( Point(point4.x-3, point4.y) )
    bridge2.append( Point(point4.x-3, point1.y) )
    bridge2.append( point1 )
    bridge2.append( point2 )
    bridge2.append( Point(point2.x+3, point2.y) )
    bridge2.append( Point(point2.x+3, point3.y) )
    bridge2.append( point3 )    
    #
    line1 = Line(Point(-100, y2), Point(100, y2))
    line2 = Line(Point(-100, y1), Point(100, y1))
    line3 = Line(Point(-100, 0.5*tw1-width7-width8-2*spacing), Point(100, 0.5*tw1-width7-width8-2*spacing)).Rotate(angle=-angle)
    line4 = Line(Point(-100, 0.5*tw1-width7-width8-width9-2*spacing), Point(100, 0.5*tw1-width7-width8-width9-2*spacing)).Rotate(angle=-angle)
    line5 = Line(Point(-100, y_5), Point(100, y_5))
    line6 = Line(Point(-100, y_6), Point(100, y_6))
    #
    point1 = line1.Intersect(line3)
    point2 = line3.Intersect(line5)
    point3 = line4.Intersect(line6)
    point4 = line2.Intersect(line4)
    #
    bridge3 = Primitive( )
    bridge3.append( point4 )
    bridge3.append( Point(point4.x-3, point4.y) )
    bridge3.append( Point(point4.x-3, point1.y) )
    bridge3.append( point1 )
    bridge3.append( point2 )
    bridge3.append( Point(point2.x+3, point2.y) )
    bridge3.append( Point(point2.x+3, point3.y) )
    bridge3.append( point3 )    
    #
    x1 = max(pt.x for pt in bridge1)
    x2 = min(pt.x for pt in bridge3)
    #
    bridge1 = bridge1.Translate(vector=-0.5*Point(x1+x2, 0))
    bridge2 = bridge2.Translate(vector=-0.5*Point(x1+x2, 0))
    bridge3 = bridge3.Translate(vector=-0.5*Point(x1+x2, 0))
    #
    return bridge1, bridge2, bridge3

#layer1 = Primitives()
#layer2 = Primitives()

#layer1.append(bridge1)
#layer2.append(bridge1.Mirror(planenormal=(1,0,0)))
#layer1.append(bridge2)
#layer2.append(bridge2.Mirror(planenormal=(1,0,0)))
#layer1.append(bridge3)
#layer2.append(bridge3.Mirror(planenormal=(1,0,0)))

