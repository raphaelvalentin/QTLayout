width = 5
spacing = 10
length = 500.0
lenght_port = 5.0

tm2 = Primitives()
poly = Primitives()


port1 = Point(lenght_port, 0.5*spacing+0.5*width, 0)
port2 = Point(lenght_port, -0.5*spacing-0.5*width, 0)
port3 = Point(lenght_port+length, 0.5*spacing+0.5*width, 0)
port4 = Point(lenght_port+length, -0.5*spacing-0.5*width, 0)



path1 = Path()
path2 = Path()

path1.append(port1)
path1.append(port3)

path2.append(port2)
path2.append(port4)



tm2.append( path1.enlarge(width=width) )
tm2.append( path2.enlarge(width=width) )

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
    





