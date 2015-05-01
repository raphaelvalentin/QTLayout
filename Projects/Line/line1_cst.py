# load the layout file
execfile('D:/Work/QTLayout/Projects/Line/line1.py')



main = CST()
main.append(Unit())
main.extend(layers.material)

mesh_metal = []
for shape in tm2:
    extrude1 = Extrude(component='tm2', name='shape1', material=metal['TM2'].material, zrange=metal['TM2'].zrange)
    extrude1.extend( shape )
    main.append( extrude1 )
    mesh_metal.append( (extrude1.component, extrude1.name) )

# port PEC
pp = Path()
pp.append( port1 )
pp.append( port1 - Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port1', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width) )
main.append( extrude1 )
mesh_metal.append( (extrude1.component, extrude1.name) )


pp = Path()
pp.append( port2 )
pp.append( port2 - Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port2', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width) )
main.append( extrude1 )
mesh_metal.append( (extrude1.component, extrude1.name) )

pp = Path()
pp.append( port3 )
pp.append( port3 + Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port3', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width) )
main.append( extrude1 )
mesh_metal.append( (extrude1.component, extrude1.name) )

pp = Path()
pp.append( port4 )
pp.append( port4 + Point(1, 0, 0) )
extrude1 = Extrude(component='port', name='port4', material='PEC', zrange=metal['TM2'].zrange)
extrude1.extend( pp.enlarge(width) )
main.append( extrude1 )
mesh_metal.append( (extrude1.component, extrude1.name) )

zport = 0.5*(metal['TM2'].zmin+metal['TM2'].zmax)
main.append( DiscretePort( P1 = port1 + Point(-1, 0, zport) ,
                       P2 = port1 + Point(-lenght_port, 0, zport) ) )
main.append( DiscretePort( P1 = port2 + Point(-1, 0, zport) ,
                       P2 = port2 + Point(-lenght_port, 0, zport) ) )
main.append( DiscretePort( P1 = port3 + Point(1, 0, zport) ,
                       P2 = port3 + Point(lenght_port, 0, zport) ) )
main.append( DiscretePort( P1 = port4 + Point(1, 0, zport) ,
                       P2 = port4 + Point(lenght_port, 0, zport) ) )


psub = Primitive()
psub.append( Point(tm2.xmin-lenght_port, tm2.ymax+50, 0) )
psub.append( Point(tm2.xmin-lenght_port, tm2.ymin-50, 0) )
psub.append( Point(tm2.xmax+lenght_port, tm2.ymin-50, 0) )
psub.append( Point(tm2.xmax+lenght_port, tm2.ymax+50, 0) )

extrude_psub = Extrude(component='env', name='psub', material=dielectric['PSUB'].material, zrange=dielectric['PSUB'].zrange)
extrude_psub.extend( psub )
main.append( extrude_psub )
extrude_air = Extrude(component='env', name='air', material=dielectric['Air'].material, zrange=dielectric['Air'].zrange)
extrude_air.extend( psub )
main.append( extrude_air )
extrude_imd1 = Extrude(component='env', name='imd1', material=dielectric['IMD1'].material, zrange=dielectric['IMD1'].zrange)
extrude_imd1.extend( psub )
main.append( extrude_imd1 )
extrude_imd2 = Extrude(component='env', name='imd2', material=dielectric['IMD2'].material, zrange=dielectric['IMD2'].zrange)
extrude_imd2.extend( psub )
main.append( extrude_imd2 )



# mesh definition
main.append( Group().Add('metal', 'mesh') )
main.append( Group().Add('imd', 'mesh') )
main.append( Group().Add('airpsub', 'mesh') )

main.append( MeshSettings('metal', 3.5)  )
main.append( MeshSettings('imd', 15)  )
main.append( MeshSettings('airpsub', 45)  )

#main.append( Group().AddItem('imd') )
main.append( Group().AddItem(extrude_psub.component, extrude_psub.name, 'airpsub') )
main.append( Group().AddItem(extrude_air.component, extrude_air.name, 'airpsub') )
main.append( Group().AddItem(extrude_imd1.component, extrude_imd1.name, 'imd') )
main.append( Group().AddItem(extrude_imd2.component, extrude_imd2.name, 'imd') )
for component, name in mesh_metal:
    main.append( Group().AddItem(component, name, 'metal') )


## Boolean Operation : Insert
for d in [extrude_imd1, extrude_imd2]:
    for m in main:
        if isinstance(m, (Brick, Solid, Extrude)) and m.component<>'env':
            op = Solid.Insert(d, m)
            main.append( op )


# simulation definitions
main.append( Boundary() )
main.append( Mesh(type='Tetrahedral') )
main.append( Solver(FrequencyRange=[0.1, 10]) )

print dielectric
