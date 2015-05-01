from syntax.__future__.path import Path
from syntax.__future__.point import Point

def Pattern(outerdiameter=360, width=5., spacing=5.):
    outerdiameter= float(outerdiameter)
    width= float(width)
    spacing= float(spacing)
    PatternShield = Primitives()
    line1 = Line(Point(-0.5*outerdiameter,0.5*outerdiameter), Point(0.5*outerdiameter,-0.5*outerdiameter))
    line1.intercept = line1.intercept + spacing
    x = 0
    for i in xrange(100):
        Block = Primitive()
        Block.append( Point(x-0.5*width,-0.5*outerdiameter) )
        Block.append( Point(x+0.5*width,-0.5*outerdiameter) )
        Block.append( Point(x+0.5*width,line1(x+width)) )
        Block.append( Point(x-0.5*width,line1(x+width)) )
        if Block.ymax-Block.ymin < width or Block.ymax > 0.5*outerdiameter:
            break
        if Block.xmax-Block.xmin < width-1e-6 or Block.xmax > 0.5*outerdiameter:
            break
        PatternShield.append( Block )
        PatternShield.append( Block.Mirror(planenormal=(0,1)) )
        if i>0:
            PatternShield.append( Block.Mirror(planenormal=(1,0)) )
            PatternShield.append( Block.Mirror(planenormal=(1,1)) )
        Block2 = Block.Rotate(center=Point(0,0), angle=90)
        PatternShield.append( Block2 )
        PatternShield.append( Block2.Mirror(planenormal=(1,0)) )
        if i>0:
            PatternShield.append( Block2.Mirror(planenormal=(0,1)) )
            PatternShield.append( Block2.Mirror(planenormal=(1,1)) )
        x = x + (width+spacing)

    Block = Primitive()
    Block.append( Point(0.5*width,0.5*outerdiameter*sqrt(2)-width) )
    Block.append( Point(-0.5*width,0.5*outerdiameter*sqrt(2)-width) )
    Block.append( Point(-0.5*width,0) )
    Block.append( Point(0.5*width,0) )
    PatternShield.append( Block.Rotate(angle=45) )
    PatternShield.append( Block.Rotate(angle=-45) )
    PatternShield.append( Block.Rotate(angle=180+45) )
    PatternShield.append( Block.Rotate(angle=180-45) )

    return PatternShield


_p = Pattern(100, 4.6, 5)