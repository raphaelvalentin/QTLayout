
outerdiameter = 150
width = 10


tm2 = Primitives()

path1 = Path()
path1.append( Point(-10, -0.5*outerdiameter) )
path1.append( Point(-0.5*outerdiameter, -0.5*outerdiameter) )
path1.append( Point(-0.5*outerdiameter, 0.5*outerdiameter) )
path1.append( Point(0, 0.5*outerdiameter) )


tm2.append( path1.enlarge(width) )
tm2.append( path1.enlarge(width).Mirror(planenormal=(1,0)) )

tm2 = tm2*1.1

