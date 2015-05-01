from exceptions import *
from math import pi, sin, cos
try:
    from CSTlib import Brick, Extrude, Point
except:
    pass
try:
    from Solid import Solid
except:
    pass
try:
    from formal import *
except:
    pass

class Rotate(object):
    def __init__(self, solid, center=(0, 0, 0), angle=(0, 0, 0)):
        if isinstance(solid, (Brick, Solid, Extrude)):
            self.solid = solid
            self.center = center
            self.angle = angle
        else:
            raise Exception('solid is not a Brick or a Solid')
    def __str__(self):
        return "\n".join([ 'With Transform',
                           '  .Reset',
                           '  .Name "{self.solid.component}:{self.solid.name}"'.format(self=self),
                           '  .Origin "Free"',
                           '  .Center "{self.center[0]}", "{self.center[1]}", "{self.center[2]}"'.format(self=self),
                           '  .Angle "{self.angle[0]}", "{self.angle[1]}", "{self.angle[2]}"'.format(self=self),
                           '  .MultipleObjects "False"',
                           '  .MultipleObjects "False"',
                           '  .GroupObjects "False"',
                           '  .Repetitions "1"',
                           '  .MultipleSelection "False"',
                           '  .Transform "Shape", "Rotate"',
                           'End With',
                           ''
                        ])
    @staticmethod
    def Rotate2D(point, center=(0, 0), angle=0):
        x, y = point
        x0, y0 = center
        theta = angle*pi/180.
        x1 = cos(theta)*(x-x0) - sin(theta)*(y-y0)+x0
        y1 = sin(theta)*(x-x0) + cos(theta)*(y-y0)+y0
        return x1, y1

    @staticmethod
    def Rotate(point, center=(0, 0, 0), angle=(0, 0, 0)):
        newPoint = Point(point.x, point.y, point.z)
        if angle[0] == angle[1] == angle[2] == 0:
            pass
        elif angle[2]<>0 and angle[1] == angle[0] == 0:
            theta = angle[2]*pi/180.
            newPoint.x = cos(theta)*(point.x-center.x) - sin(theta)*(point.y-center.y)+center.x
            newPoint.y = sin(theta)*(point.x-center.x) + cos(theta)*(point.y-center.y)+center.y
        elif angle[0]<>0 and angle[1] == angle[2] == 0:
            theta = angle[0]*pi/180.
            raise Exception('Rotate operation not implemented.')
        elif angle[1]<>0 and angle[0] == angle[2] == 0:
            theta = angle[1]*pi/180.
            raise Exception('Rotate operation not implemented.')
        else:
            raise Exception('Rotate around one angle a time.')
        return newPoint
    
    def __copy__(self):
        return Rotate( solid=self.solid, center=self.center, angle=self.angle)


class Mirror(object):
    def __init__(self, solid, center=(0, 0, 0), planenormal=(0, 0, 0)):
        if isinstance(solid, (Brick, Solid, Extrude)):
            self.solid = solid
            self.center = center
            self.planenormal = planenormal
        else:
            raise Exception('solid is not a Brick or a Solid')
    def __str__(self):
        return "\n".join([  'With Transform',
                          '  .Reset',
                          '  .Name "{self.solid.component}:{self.solid.name}"'.format(self=self),
                          '  .Origin "Free"',
                          '  .Center "{self.center[0]}", "{self.center[1]}", "{self.center[2]}"'.format(self=self),
                          '  .PlaneNormal "{self.planenormal[0]}", "{self.planenormal[1]}", "{self.planenormal[2]}"'.format(self=self),
                          '  .MultipleObjects "False"',
                          '  .GroupObjects "False"',
                          '  .Repetitions "1"',
                          '  .MultipleSelection "False"',
                          '  .Destination ""',
                          '  .Material ""',
                          '  .Transform "Shape", "Mirror"',
                          '  End With',
                          ''
                       ])

class Translate(object):
    def __init__(self, solid, vector=(0, 0, 0)):
        if isinstance(solid, (Brick, Solid, Extrude)):
            self.solid = solid
            self.vector = vector
        else:
            raise Exception('solid is not a Brick or a Solid')
    def __str__(self):
        return "\n".join([ 'With Transform',
                           '  .Reset',
                           '  .Name "{self.solid.component}:{self.solid.name}"'.format(self=self),
                           '  .Vector "{self.vector[0]}", "{self.vector[1]}", "{self.vector[2]}"'.format(self=self),
                           '  .UsePickedPoints "False" ',
                           '  .InvertPickedPoints "False" ',
                           '  .MultipleObjects "False"',
                           '  .GroupObjects "False"',
                           '  .Repetitions "1"',
                           '  .MultipleSelection "False" ',
                           '  .Transform "Shape", "Translate"',
                           '  End With',
                           ''
                       ])

