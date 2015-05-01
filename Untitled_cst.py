# load the layout file
execfile('.\Projects\Transformer-PA-2nH-Ratio2_4_Edward_run2/transformer10.py')
from syntax.smic import *

main = CST() 
main.append(Unit())
main.extend(material)
for shape in tm2:
    extrude1 = Extrude(component='tm2', name='shape1', material=metal['TM2'].material, zrange=metal['TM2'].zrange)
    extrude1.extend( shape )
    main.append( extrude1 )
for shape in alrdl:
    extrude1 = Extrude(component='alrdl', name='shape2', material=metal['ALRDL'].material, zrange=metal['ALRDL'].zrange)
    extrude1.extend( shape )
    main.append( extrude1 )
zrange = min(metal['Metal1'].zrange), max(metal['Metal4'].zrange)
for shape in metal1:
    extrude1 = Extrude(component='metal1234', name='shape2', material=metal['Metal1'].material, zrange=zrange)
    extrude1.extend( shape )
    main.append( extrude1 )




inter1 = Layer(name='Inter1', material='Aluminum', zmin=min(metal['ALRDL'].zrange), zmax=max(metal['TM2'].zrange))
inter2 = Layer(name='Inter2', material=metal['Metal6'].material, zmin=min(metal['TM2'].zrange), zmax=max(metal['Metal6'].zrange))
inter3 = Layer(name='Inter3', material=metal['Metal5'].material, zmin=min(metal['Metal6'].zrange), zmax=max(metal['Metal5'].zrange))
inter4 = Layer(name='Inter4', material=metal['Metal4'].material, zmin=min(metal['Metal5'].zrange), zmax=max(metal['Metal4'].zrange))


for shape in contact_tm2_alrdl:
    extrude1 = Extrude(component='via1', name='shape1', material=metal['ALRDL'].material, zrange=(metal['ALRDL'].zmin, metal['TM2'].zmax) )
    extrude1.extend( shape )
    main.append( extrude1 )



try: 
    for shape in poly:
        extrude1 = Extrude(component='metal1', name='shape1', material='PEC', zrange=metal['Poly'].zrange)
        extrude1.extend( shape )
        main.append( extrude1 )
except:
    pass

# port
zrange = metal['TM2'].zmin, metal['ALRDL'].zmax
extrude1 = Extrude(component='port', name='port1', material='PEC', zrange=zrange)
extrude1.append( port1[0] )
extrude1.append( port1[1] )
extrude1.append( port1[1]+Point(0, 1) )
extrude1.append( port1[0]+Point(0, 1) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port2', material='PEC', zrange=zrange)
extrude1.append( port2[0] )
extrude1.append( port2[1] )
extrude1.append( port2[1]+Point(0, 1) )
extrude1.append( port2[0]+Point(0, 1) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port3', material='PEC', zrange=zrange)
extrude1.append( port3[0] )
extrude1.append( port3[1] )
extrude1.append( port3[1]-Point(0, 1) )
extrude1.append( port3[0]-Point(0, 1) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port4', material='PEC', zrange=zrange)
extrude1.append( port4[0] )
extrude1.append( port4[1] )
extrude1.append( port4[1]-Point(0, 1) )
extrude1.append( port4[0]-Point(0, 1) )
main.append( extrude1 )

# env
extrude1 = Extrude(component='env', name='si', material=dielectric['PSUB'].material, zrange=dielectric['PSUB'].zrange)
extrude1.append( Point(poly.xmin-30, poly.ymin-30) )
extrude1.append( Point(poly.xmax+30, poly.ymin-30) )
extrude1.append( Point(poly.xmax+30, poly.ymax+30) )
extrude1.append( Point(poly.xmin-30, poly.ymax+30) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port1', material='PEC', zrange=zrange)
extrude1.append( Point(port1[0].x, poly.ymax+30) )
extrude1.append( Point(port1[1].x, poly.ymax+30) )
extrude1.append( port1[1]+Point(0, 5) )
extrude1.append( port1[0]+Point(0, 5) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port2', material='PEC', zrange=zrange)
extrude1.append( Point(port2[0].x, poly.ymax+30) )
extrude1.append( Point(port2[1].x, poly.ymax+30) )
extrude1.append( port2[1]+Point(0, 5) )
extrude1.append( port2[0]+Point(0, 5) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port3', material='PEC', zrange=zrange)
extrude1.append( Point(port3[0].x, poly.ymin-30) )
extrude1.append( Point(port3[1].x, poly.ymin-30) )
extrude1.append( port3[1]-Point(0, 5) )
extrude1.append( port3[0]-Point(0, 5) )
main.append( extrude1 )

extrude1 = Extrude(component='port', name='port4', material='PEC', zrange=zrange)
extrude1.append( Point(port4[0].x, poly.ymin-30) )
extrude1.append( Point(port4[1].x, poly.ymin-30) )
extrude1.append( port4[1]-Point(0, 5) )
extrude1.append( port4[0]-Point(0, 5) )
main.append( extrude1 )



# simulation definitions
main.append( Boundary() )
main.append( Mesh(type='Tetrahedral') )
main.append( Solver(FrequencyRange=[0.1, 10]) )

print dielectric

