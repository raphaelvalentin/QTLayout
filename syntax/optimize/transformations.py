from math import cos, acos, pi, isinf, isnan, sqrt

__all__ = ['BoxContraintsTransformation', 'ScalingTransformation']

nan = float('nan')
inf = float('inf')
       
class BoxContraintsTransformation:
    def __init__(self, bounds=()):
        bounds = tuple(bounds)
	 
        if bounds == ():
            xmin, xmax = nan, nan
	else:
            xmin, xmax = float(bounds[0]), float(bounds[1])
	    
        if xmax<xmin:
           xmin, xmax = xmax, xmin
	   
        if not isinf(xmin) and not isinf(xmax) and not isnan(xmin) and not isnan(xmax) and xmin==xmax:
            def f(x): 
                return xmin
            def finv(y): 
                return xmin

        elif ( isnan(xmin) and isnan(xmax) ) or ( xmin==-inf and xmax==inf ):
            def f(x):
                return x
            def finv(y): 
                return y
    	    
        elif not isinf(xmin) and not isinf(xmax) and not isnan(xmin) and not isnan(xmax):
            def f(x): 
                return xmin + 0.5*(xmax-xmin) * ( 1.0 - cos( (x-xmin)*pi/(xmax-xmin) ) )
            def finv(y): 
                return (xmax-xmin)/pi * acos( 1.0 - 2.0*(y-xmin)/(xmax-xmin) ) + xmin
    	    
        elif not isinf(xmin) and (xmax==inf or isnan(xmax)):
            def f(x): 
                return xmin + (x-xmin)**2
            def finv(y): 
                return sqrt(y-xmin) + xmin
            
        elif not isinf(xmax) and (xmin==-inf  or isnan(xmin)):
            def f(x): 
                return xmax - (x-xmax)**2
            def finv(y): 
                return -sqrt(xmax-y) + xmax
    
        else:
            raise Exception('Boundary conditions are not set up correctly.')

        self.__call__ = f
	self.inverse = finv
        

class ScalingTransformation:
    def __init__(self, scaling):
         scaling = float(scaling)

	 if scaling>0:
   	     def f(x):
	         return x/scaling
	     def finv(y):
	         return y*scaling
         else:
            raise Exception('Scaling Transformation is not set up correctly.')

         self.__call__ = f
	 self.inverse = finv
	 
        	 

    
    
