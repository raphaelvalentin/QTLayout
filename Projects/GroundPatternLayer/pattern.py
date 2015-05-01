
def Pattern(outerdiameter=360, width=5., spacing=5.):
    PatternShield = Primitives()
    line1 = Line(Point(-0.5*outerdiameter,0.5*outerdiameter), Point(0.5*outerdiameter,-0.5*outerdiameter))
    v = Point(0, spacing)
    #line1.intercept = line1.intercept + spacing
    line1 = line1.Translate(vector=v)
    x = 0
    for i in xrange(100):
        if line1(x+0.5*width)+0.5*outerdiameter < width:
            break
        Block = Primitive()
        Block.append( Point(x-0.5*width,-0.5*outerdiameter) )
        Block.append( Point(x+0.5*width,-0.5*outerdiameter) )
        Block.append( Point(x+0.5*width,line1(x+0.5*width)) )
        Block.append( Point(x-0.5*width,line1(x+0.5*width)) )
        PatternShield.append( Block )
        PatternShield.append( Block.Mirror(planenormal=(1,0)) )
        PatternShield.append( Block.Mirror(planenormal=(0,1)) )
        PatternShield.append( Block.Mirror(planenormal=(1,1)) )
        Block2 = Block.Rotate(center=Point(0,0), angle=90)
        PatternShield.append( Block2 )
        PatternShield.append( Block2.Mirror(planenormal=(1,0)) )
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

#poly = Pattern( outerdiameter=370+2*30 )

class Polygon(list):
    def __init__(self, diameter, num):
        self.num = num
        self.diameter = diameter
    def __getitem__(self, i):
        theta = pi/self.num
        rayon = 0.5*self.diameter/cos(theta)
        if i==0:
            return Point( 0.5*self.diameter, 0.5*self.diameter*tan(theta), 0.0 )
        else:
            return rayon*Point( cos((2*i+1)*theta), sin((2*i+1)*theta), 0.0 )


def Pattern2(outerdiameter=360, width=5., spacing=5.):
    PatternShield = Primitives()
    Polygon1 = Primitive()
    Polygon1.append( Point(0.5*outerdiameter, 0, 0) )
    for i in xrange(12/4):
        Polygon1.append( Polygon(outerdiameter, 12)[i] )
    Polygon1.append( Point(0, 0.5*outerdiameter,  0) )
    Polygon1 = Polygon1.Rotate(angle=(0,0,-90))
    line1 = Line(Point(-0.5*outerdiameter,0.5*outerdiameter,0), Point(0.5*outerdiameter,-0.5*outerdiameter,0))
    line1.intercept = line1.intercept + spacing
    x = 0
    for i in xrange(100):
        if line1(x+0.5*width)+0.5*outerdiameter < width:
            break
        Block = Primitive()
        ymin = -0.5*outerdiameter
        for j in xrange(len(Polygon1)):
            if Polygon1[j].x<=x<=Polygon1[j+1].x:
                line2 = Line( Polygon1[j], Polygon1[j+1] )
                ymin = line2.slope * x + line2.intercept
                break
        if ymin-line1(x+0.5*width)>=0:
            break
        Block.append( Point(x-0.5*width,ymin,0) )
        Block.append( Point(x+0.5*width,ymin,0) )
        Block.append( Point(x+0.5*width,line1(x+0.5*width),0) )
        Block.append( Point(x-0.5*width,line1(x+0.5*width),0) )
        PatternShield.append( Block )
        PatternShield.append( Block.Mirror(planenormal=(1,0,0)) )
        PatternShield.append( Block.Mirror(planenormal=(0,1,0)) )
        PatternShield.append( Block.Mirror(planenormal=(1,1,0)) )
        Block2 = Block.Rotate(center=Point(0,0,0), angle=(0,0,90))
        PatternShield.append( Block2 )
        PatternShield.append( Block2.Mirror(planenormal=(1,0,0)) )
        PatternShield.append( Block2.Mirror(planenormal=(0,1,0)) )
        PatternShield.append( Block2.Mirror(planenormal=(1,1,0)) )
        x = x + (width+spacing)

    
    Path1 = Path()
    Path1.append( Point(0, 0, 0) )
    Path1.append( Point(Polygon1[2].x+4, Polygon1[2].y-4, 0) )
    Block = Path1.enlarge(width=width)
    PatternShield.append( Block )
    PatternShield.append( Block.Rotate(angle=(0,0,-90)) )
    PatternShield.append( Block.Rotate(angle=(0,0,180)) )
    PatternShield.append( Block.Rotate(angle=(0,0,90)) )

    Polygon2 = Primitive()
    Polygon3 = Primitive()
    Polygon2.append( Point(0.5*outerdiameter+3, 0.5*(spacing+width), 0) )
    Polygon3.append( Point(0.5*outerdiameter+3+width, 0.5*(spacing+width), 0) )
    for i in xrange(12/4):
        Polygon2.append( Polygon(outerdiameter+6, 12)[i] )
        Polygon3.append( Polygon(outerdiameter+6+2*width, 12)[i] )
    Polygon2.append( Point(0.5*(spacing+width), 0.5*outerdiameter+3,  0) )
    Polygon3.append( Point(0.5*(spacing+width), 0.5*outerdiameter+3+width,  0) )

    Ring = Primitive()
    Ring.extend( Polygon2 )
    Ring.extend( list(reversed(Polygon3)) )

    PatternShield.append( Ring )
    PatternShield.append( Ring.Rotate(angle=(0,0,-90)) )
    PatternShield.append( Ring.Rotate(angle=(0,0,180)) )
    PatternShield.append( Ring.Rotate(angle=(0,0,90)) )



    return PatternShield

#main = CST()
#pattern = Primitives()
#pattern.extend( Pattern(outerdiameter=280+40).Translate(vector=Point(0.0, 121.24355653, 0.0)) )
#pattern.extend( Pattern(outerdiameter=280+40).Translate(vector=Point(0.0, -121.24355653, 0.0)) )
#main.extend(pattern)