width = 16/1.1
spacing = 9.5/1.1
length = 500.0
lenght_port = 5.0

tm2 = Primitives()
poly = Primitives()

spacing2 = 3.0

port1 = Point(lenght_port, 0.5*spacing+0.5*width+0.5*spacing2, 0)
port2 = Point(lenght_port, -0.5*spacing-0.5*width-0.5*spacing2, 0)
port3 = Point(lenght_port+length, 0.5*spacing+0.5*width+0.5*spacing2, 0)
port4 = Point(lenght_port+length, -0.5*spacing-0.5*width-0.5*spacing2, 0)



path1a = Path()
path1b = Path()
path2a = Path()
path2b = Path()

path1a.append(port1+Point(0,0.5*spacing2+0.25*width,0))
path1a.append(port3+Point(0,0.5*spacing2+0.25*width,0))
path1b.append(port1-Point(0,0.5*spacing2+0.25*width,0))
path1b.append(port3-Point(0,0.5*spacing2+0.25*width,0))

path2a.append(port2+Point(0,0.5*spacing2+0.25*width,0))
path2a.append(port4+Point(0,0.5*spacing2+0.25*width,0))
path2b.append(port2-Point(0,0.5*spacing2+0.25*width,0))
path2b.append(port4-Point(0,0.5*spacing2+0.25*width,0))



tm2.append( path1a.enlarge(width=0.5*width) )
tm2.append( path1b.enlarge(width=0.5*width) )
tm2.append( path2a.enlarge(width=0.5*width) )
tm2.append( path2b.enlarge(width=0.5*width) )

width_poly = 5.0
spacing_poly = 5.0
length_poly = width*2+spacing + 4*spacing

x = 5.0+0.5*width_poly
for i in xrange(1000):
    if x>length:
        break
    path = Path()
    path.append( Point(x, 0.5*length_poly, 0  ) )
    path.append( Point(x, -0.5*length_poly, 0  ) )
    poly.append( path.enlarge(width_poly) )
    x += width_poly+spacing_poly
    





