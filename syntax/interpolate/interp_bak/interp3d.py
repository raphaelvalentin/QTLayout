from libarray import *
from exceptions import Exception
from linalg import solve, det

def interp3d(x, y, z, data):

    class point(list):
        def __init__(self, *args):
	    list.__init__(self, args)
    pt0 = point(x, y, z)

    # check that the ref point is inside the data points
    xs, ys, zs, fs = data.T
    if not(min(xs)<=x<=max(xs)) or not(min(ys)<=y<=max(ys)) or not(min(zs)<=z<=max(zs)):
        raise Exception('Extrapolation is not allowed.')
    
    def dist(pt0, pt1):
        # compute the absolute distance
	try:
            return sqrt(sum([ (2*(a-b)/(a+b))**2 for a, b in zip(pt0, pt1)]))
	except:
            return sqrt(sum([ (a-b)**2 for a, b in zip(pt0, pt1) ]))
    # generate a list of distance between the ref point and each data points
    dist0 = [ (i, dist(point(x, y, z), pt0)) for i, (x, y, z, f) in enumerate(data) ]
    # sort
    dist0 = sorted(dist0, key=lambda x : x[1])
    
    if dist0[1]==0.0:
	return data[dist0[0]][-1]
	
    # initialise
    A = array([])
    B = array([])

    # first point
    x, y, z, f = data[dist0[0][0]]
    # build the system a*x+b*y+c*z+d = f -> linear interpolation
    A.append([x, y, z, 1.])
    B.append([f])
    
    # second point
    for k, v in dist0:
        x, y, z, f = data[k]
	if not(x in A.T[0]):
            A.append([x, y, z, 1.])
            B.append([f])
	    break
	    
    # third point
    for k, v in dist0:
        x, y, z, f = data[k]
	if not(y in A.T[1]):
            A.append([x, y, z, 1.])
            B.append([f])
	    break
	    
    # fourth point
    for k, v in dist0:
        x, y, z, f = data[k]
	if not(z in A.T[2]):
            A.append([x, y, z, 1.])
            B.append([f])
	    break
	    
    X = solve(A, B)
    a, b, c, d, = X.T[0]
    
    # check that the system give accurate result
    err = 0.0; ref = 0.0
    for [x, y, z, dummy], f in zip(A, B.T[0]):
        err += (a*x+b*y+c*z+d - f)**2
	ref += f*f
    if ref<>0:
        err = sqrt(err/ref)
    else:
        err = sqrt(err)
    if err>1e-6:
        Exception('The linear interpolation is not accurate.')
    
    return a*pt0[0]+b*pt0[1]+c*pt0[2]+d

if __name__ == '__main__':
    
    """"
    from spectre.lib.smic.smic40ll import *
    filename='//home/vraphael/.python/spectre/lib/smic/smic40ll/layout_effect.dat'
    from rawdata import *
    layout_effect = base(filename).read()
    t = layout_effect[0]
    core = array([t['lr'], t['wr'], t['nf'], t['scar']]).T
    print interp3d(0.05e-6, 2.0e-6, 32, core)
    
    #print n11ll_ckt_rf(l=0.04e-6, w=2e-6, nf=32)['scar']
    #print n11ll_ckt_rf(l=0.041e-6, w=2e-6, nf=32)['scar']
    exit(0)
    """

    x = array(range(-10,10,1), dtype=float)
    y = array(range(-10,10,1), dtype=float)
    z = array(range(-10,10,1), dtype=float)

    data = array([])


    def f(x, y, z):
        return ((x)**2+y**2+z**2)
    
    for xi in x:
        for yi in y:
            for zi in z:
	        data.append([xi, yi, zi, f(xi, yi, zi)])


    print interp3d(1.1, 1.5, 2.5, data)
    print f(1.1, 1.5, 2.5)
     
