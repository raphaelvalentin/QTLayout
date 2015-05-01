# load the layout file
execfile('D:/Work/QTLayout/Projects/Line/line2.py')



main = CST()
main.append(Unit())
main.extend(layers.material)
for shape in tm2:
    extrude1 = Extrude(component='tm2', name='shape1', material=metal['TM2'].material, zrange=metal['TM2'].zrange)
    extrude1.extend( shape )
    main.append( extrude1 )
for shape in poly:
    extrude1 = Extrude(component='poly', name='shape1', material=metal['Poly'].material, zrange=metal['Poly'].zrange)
    extrude1.extend( shape )
    main.append( extrude1 )

# port PEC
pp = Path()
pp.append( port1 )
pp.append( port1 - Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port1', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width+spacing2) )
main.append( extrude1 )

pp = Path()
pp.append( port2 )
pp.append( port2 - Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port2', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width+spacing2) )
main.append( extrude1 )

pp = Path()
pp.append( port3 )
pp.append( port3 + Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port3', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width+spacing2) )
main.append( extrude1 )

pp = Path()
pp.append( port4 )
pp.append( port4 + Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port4', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width+spacing2) )
main.append( extrude1 )

zport = 0.5*(metal['TM2'].zmin+metal['TM2'].zmax)
port1_ = DiscretePort( P1 = port1 + Point(-1, 0, zport) ,
                       P2 = port1 + Point(-lenght_port, 0, zport) )
port2_ = DiscretePort( P1 = port2 + Point(-1, 0, zport) ,
                       P2 = port2 + Point(-lenght_port, 0, zport) )
port3_ = DiscretePort( P1 = port3 + Point(1, 0, zport) ,
                       P2 = port3 + Point(lenght_port, 0, zport) )
port4_ = DiscretePort( P1 = port4 + Point(1, 0, zport) ,
                       P2 = port4 + Point(lenght_port, 0, zport) )

main.append( port1_ )
main.append( port2_ )
main.append( port3_ )
main.append( port4_ )

psub = Primitive()
psub.append( Point(tm2.xmin-lenght_port, tm2.ymax+50, 0) )
psub.append( Point(tm2.xmin-lenght_port, tm2.ymin-50, 0) )
psub.append( Point(tm2.xmax+lenght_port, tm2.ymin-50, 0) )
psub.append( Point(tm2.xmax+lenght_port, tm2.ymax+50, 0) )

extrude1 = Extrude(component='env', name='psub', material=dielectric['PSUB'].material, zrange=dielectric['PSUB'].zrange)
extrude1.extend( psub )
main.append( extrude1 )
extrude1 = Extrude(component='env', name='air', material=dielectric['Air'].material, zrange=dielectric['Air'].zrange)
extrude1.extend( psub )
main.append( extrude1 )


# simulation definitions
main.append( Boundary() )
main.append( Mesh(type='Tetrahedral') )
main.append( Solver(FrequencyRange=[0.1, 10]) )

print dielectric
