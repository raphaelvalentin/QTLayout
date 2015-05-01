import gdspy


class Cell(gdspy.Cell):
    layer = 1
    def __init__(self, name):
        gdspy.Cell.__init__(self, name)
    def append(self, primitives, layer=None):
        if layer:
            self.layer = layer
        else:
            self.layer = Cell.layer
        for primitive in primitives:
            points = [(pt.x, pt.y) for pt in primitive]
            poly1 = gdspy.Polygon(points, layer)
            self.add(poly1)
            
        
global __unit__
__unit__ = 1.0e-6
def unit(x):
    global __unit__
    __unit__=x
    
global __precision__
__precision__ = 5.0e-9
def precision(x):
    global __precision__
    __precision__=x
    


def export(filename):
    global __precision__
    global __unit__
    gdspy.gds_print(filename, unit=__unit__, precision=__precision__)

def viewer():
    gdspy.LayoutViewer()

        
