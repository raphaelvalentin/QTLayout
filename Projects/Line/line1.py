width = 16/1.1
spacing = 9.5/1.1

width = 5.0
spacing = 5.0
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





