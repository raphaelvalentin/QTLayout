from libarray import *
from function import trunc
from LinearNDInterpolator import LinearNDInterpolator
from Quadratic1DInterpolate_C1 import Quadratic1DInterpolate
from exceptions import Exception

__all__ = ['LinearNDInterpolator', 'lookup', 'Quadratic1DInterpolate']

def lookup(x, y, z, table, mod='column'):
    if mod == 'column':
        table = table.T
    points, mvalues = array([table[0], table[1], table[2]]).T, table[3:]
    result = array([x, y, z])
    for values in mvalues:
        try:
	    LI = LinearNDInterpolator(points, values)
     	    result.append(trunc(LI((x, y, z)), 14))
	except:
	    raise Exception('Error in lookup function with x, y, z = %g, %g, %g and values = %r'%(x, y, z, values))
    return result


if __name__ == '__main__':
    filename='/home/vraphael/.python/spectre/lib/tsmc/crn40lp/layout_effect.dat'
    from rawdata import *
    layout_effect = base(filename).read()
    t = layout_effect[0]
    core = array([t['lr'], t['wr'], t['nr'], t['spba'], t['spa2'], t['spa3'], t['spa1'], t['sca'], t['scc'], t['scb'], t['spa'], t['sap'], t['sodx1'], t['sodx2'], t['spba1'], t['sa3'], t['sa2'], t['sa1'], t['sodx'], t['sa6'], t['sa5'], t['sa4'], t['sapb'], t['sody'], t['sa']])
    print lookup(0.041e-6, 2e-6, 16, core, mod='row')

