class ClosedCurveError(Exception):
    pass


class Point(object):
    def __init__(self, x, y, z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    def __str__(self):
        return "(%r %r %r)"%(self.x, self.y, self.z)
    def __repr__(self):
        return "(%r %r %r)"%(self.x, self.y, self.z)
    def __eq__(self, point):
        return self.x==point.x and self.y==point.y and self.z == point.z

class ThickenSheet(object):
    unit = 'um'
    def __init__(self, **kwargs):
        self.sheet = kwargs.get('sheet', 'sheet1')
        self.thickness = kwargs.get('thickness', 1.0)
    def __str__(self):
        s = ['oEditor.ThickenSheet',
            '  Array("NAME:Selections",',
            '    "Selections:=", "{sheet}",'.format(sheet=self.sheet),
            '    "NewPartsModelFlag:=", "Model" ),',
            '  Array("NAME:SheetThickenParameters",',
            '    "Thickness:=", "{thickness}{unit}",'.format(thickness=self.thickness, unit=ThickenSheet.unit),
            '    "BothSides:=", false )\n']
        return " _\n".join(s)

class CreatePolyline(object):
    unit = 'um'
    def __init__(self, *args, **kwargs):
        self.points = args
        self.name = kwargs.get('name', 'Polyline1')
        self.material = kwargs.get('material', 'vacuum')
        self.check()
    
    @property
    def edges(self):
        for i in xrange(len(self.points)-1):
            yield (self.points[i], self.points[i+1])
        yield (self.points[i+1], self.points[0])
    
    def check(self):
        for i, section1 in enumerate(self.edges):
            for j, section2 in enumerate(self.edges):
                if i<>j:
                    if section1[0]==section2[0] and section1[1]==section2[1] or \
                       section1[1]==section2[0] and section1[0]==section2[1]:
                        raise ClosedCurveError('there is closed curve in the list.')

    def __str__(self):
        s = []
        s.append( 'oEditor.CreatePolyline' )
        s.append( '  Array("NAME:PolylineParameters",' )
        s.append( '    "IsPolylineCovered:=", true,' )
        s.append( '    "IsPolylineClosed:=", true,' )
        s.append( '    Array("NAME:PolylinePoints",' )
        for point in self.points:
            s.append( '      Array("NAME:PLPoint", "X:=", "{x}{unit}", "Y:=", "{y}{unit}", "Z:=", "{z}{unit}"),'.format( x=point.x, y=point.y, z=point.z, unit=CreatePolyline.unit ) )
        s.append( '      Array("NAME:PLPoint", "X:=", "{x}{unit}", "Y:=", "{y}{unit}", "Z:=", "{z}{unit}")),'.format( x=self.points[0].x, y=self.points[0].y, z=self.points[0].z, unit=CreatePolyline.unit ) )
        
        s.append( '    Array("NAME:PolylineSegments",' )
        for i in xrange(len(self.points)-1):
            s.append( '      Array("NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", {indx}, "NoOfPoints:=", 2),'.format( indx=i ) )
        s.append( '      Array("NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", {indx}, "NoOfPoints:=", 2)),'.format( indx=len(self.points)-1 ) )
        s.append( '    Array(' )
        s.append( '      "NAME:PolylineXSection",' )
        s.append( '      "XSectionType:=", "None",' )
        s.append( '      "XSectionOrient:=", "Auto",' )
        s.append( '      "XSectionWidth:=", "0{unit}",'.format( unit=CreatePolyline.unit ) )
        s.append( '      "XSectionTopWidth:=", "0{unit}",'.format( unit=CreatePolyline.unit ) )
        s.append( '      "XSectionHeight:=", "0{unit}",'.format( unit=CreatePolyline.unit ) )
        s.append( '      "XSectionNumSegments:=", "0",' )
        s.append( '      "XSectionBendType:=", "Corner" )),' )
        s.append( '  Array("NAME:Attributes",'  )
        s.append( '    "Name:=", "{name}",'.format( name=self.name ) )
        s.append( '    "Flags:=", "",' )
        s.append( '    "Color:=", "(132 132 193)",' )
        s.append( '    "Transparency:=", 0,' )
        s.append( '    "PartCoordinateSystem:=", "Global",' )
        s.append( '    "UDMId:=", "",' )
        s.append( '    "MaterialValue:=", "" & Chr(34) & "{material}" & Chr(34) & "",'.format( material = self.material) )
        s.append( '    "SolveInside:=", true )\n' )
        return " _\n".join(s)
        
class Extrude(list):
    unit = 'um'
    def __init__(self, **kwargs):
        self.name = kwargs.get('name', 'Extrude')
        self.component = kwargs.get('component', 'component1')
        self.material = kwargs.get('material', 'vacuum')
        self.zrange = kwargs.get('zrange', (0.0, 1.0))
    def __str__(self):
        s = []
        points = (Point(point.x, point.y, self.zrange[0]) for point in self)
        s.append( str(CreatePolyline(*points, name=self.name, material=self.material)) )
        s.append( str(ThickenSheet(sheet=self.name, thickness=self.zrange[0]-self.zrange[1])) )
        return "".join(s)

class Material:
    tnom = 25.0
    temp = 25.0
    def __init__(self, name='material1', **kwargs):
        self.name = "{name}".format(name=name)
        self.epsilon = kwargs.get('epsilon', 1.0)
        self.conductivity = kwargs.get('kappa', 0.0)
        self.type  = kwargs.get('type', "Normal")
        self.tc1 = kwargs.get('tc1', 0.0)
    def __str__(self):
        s.append( 'oDefinitionManager.AddMaterial' )
        s.append( '  Array("NAME:{name}",'.formal(name=self.name) )
        s.append( '    "CoordinateSystemType:=", "Cartesian",' )
        s.append( '    Array("NAME:AttachedData"),' )
        s.append( '    Array("NAME:ModifierData"),' )
        s.append( '    "conductivity:=", "{conductivity}")'.format(conductivity=self.conductivity) )
        return " _\n".join(s)


class HFSS(list):
    Formatting = ":.5f"
    def __init__(self, **kwargs):
        self.component = kwargs.get('component', 'component1')
        self.name = kwargs.get('name', 'object1')
    def __str__(self):
        s = []
        s.append( 'Dim oAnsoftApp ')
        s.append( 'Dim oDesktop ')
        s.append( 'Dim oProject ')
        s.append( 'Dim oDesign ')
        s.append( 'Dim oEditor ')
        s.append( 'Dim oModule ')
        s.append( 'Set oAnsoftApp = CreateObject("AnsoftHfss.HfssScriptInterface") ')
        s.append( 'Set oDesktop = oAnsoftApp.GetAppDesktop() ')
        s.append( 'oDesktop.RestoreWindow ')
        s.append( 'Set oProject = oDesktop.GetActiveProject() ')
        s.append( 'Set oDefinitionManager = oProject.GetDefinitionManager() ')
        s.append( 'Set oDesign = oProject.SetActiveDesign("HFSSDesign1") ')
        s.append( 'Set oEditor = oDesign.SetActiveEditor("3D Modeler") ')
        for line in self:
            if not isinstance(line, HFSS):
                s.append( str(line) )
            else:
                s.append( "\n".join(str(line).split('\n')[2:-1]) )
        return "\n".join(s)


a = []
a.append( Point(2, 0.1) )
a.append( Point(2, -2) )
a.append( Point(-2, -2) )
a.append( Point(-2, 2) )
a.append( Point(2, 2) )
a.append( Point(2, 0) )
a.append( Point(1, 0) )
a.append( Point(1, 1) )
a.append( Point(-1, 1) )
a.append( Point(-1, -1) )
a.append( Point(1, -1) )
a.append( Point(1, 0.1) )

b = Extrude(zrange=(0.5, 0.6))
b.extend(a)

project1 = HFSS()
project1.append(b)

print project1
