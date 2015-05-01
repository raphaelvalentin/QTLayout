import matplotlib.pyplot as plt
from CSTlib import *
from newton import NewtonRaphson as fmin
from libind import *


######################################
# parameters loops
width = 10
outerdiameter = 250.0
n = 112

spacing = 7.0

Tnom = 25.
T = 25.
tc1 = 0.0039


layer1 = Primitives()
circle1, circle2 = Ring(outerdiameter=outerdiameter, width=width , n=n)

segment1 = Segment( Point(-0.5*spacing,-outerdiameter,0), Point(-0.5*spacing,0,0) )

# first section
block2 = Primitive()
for segment in circle1.edges:
    if not segment.Intercept(segment1):
        block2.append( segment[0] )
    else:
        point1 = Line(*segment1).Intercept(Line(*segment))
        block2.append( point1 )
        break

isIntercept = False
for segment in reversed(circle2.edges):
    if segment.Intercept(segment1):
        isIntercept = True
        point1 = Line(*segment1).Intercept(Line(*segment))
        block2.append( point1 )
    elif isIntercept:
        block2.append( segment[0] )

# second section        
segment2 = Segment( Point(0.5*spacing,-outerdiameter,0), Point(0.5*spacing,0,0) )
block3 = Primitive()
for segment in reversed(circle1.edges):
    if not segment.Intercept(segment2):
        block3.append( segment[1] )
    else:
        point1 = Line(*segment2).Intercept(Line(*segment))
        block3.append( point1 )
        break

isIntercept = False
for segment in circle2.edges:
    if segment.Intercept(segment2):
        isIntercept = True
        point1 = Line(*segment2).Intercept(Line(*segment))
        block3.append( point1 )
    elif isIntercept:
        block3.append( segment[1] )



layer1.append( block2 )
layer1.append( block3 )

for shape in layer1:
    plt.plot(list(point.x for point in shape), list(point.y for point in shape), color='black')
    plt.fill(list(point.x for point in shape), list(point.y for point in shape), color='red')
    
plt.axes().set_aspect('equal', 'datalim')
plt.show()



from layers import *
import layers
main = CST()
main.append(Unit())
main.extend(layers.material)
material.append( Material(name='Aluminum', kappa=3.831418e+07/(1.0+tc1*(T-Tnom)), color=(0.75, 0.75, 0.75)) )
tm2 = metal['TM2']
alrdl = metal['ALRDL']

for shape in layer1:
    extrude1 = Extrude(component='component2', name='alrdl', material=alrdl.material, zrange=alrdl.zrange)
    extrude1.extend( shape )
    main.append( extrude1 )

    
    
main.append(Boundary())
main.append(Mesh() )
main.append(Solver(FrequencyRange=(0.1,50)) )

main.append("""Group.Add "metal", "mesh"
                """)
main.append("""Group.Add "imd", "mesh"
                """)
main.append("""Group.Add "airsi", "mesh"
                """)


main.append("""With MeshSettings
     With .ItemMeshSettings ("group$metal")
          .SetMeshType "Tet"
          .Set "OctreeSizeFaces", "0"
          .Set "Size", "2"
     End With
End With
""")
    
main.append("""With MeshSettings
     With .ItemMeshSettings ("group$imd")
          .SetMeshType "Tet"
          .Set "OctreeSizeFaces", "0"
          .Set "Size", "10"
     End With
End With
""")
    
main.append("""With MeshSettings
     With .ItemMeshSettings ("group$airsi")
          .SetMeshType "Tet"
          .Set "OctreeSizeFaces", "0"
          .Set "Size", "45"
     End With
End With
""")
  
        
print main




