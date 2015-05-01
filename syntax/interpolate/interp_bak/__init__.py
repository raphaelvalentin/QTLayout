from libarray import *
from interp2d import interp2d
from interp3d import interp3d
from function import trunc

__all__ = ['interp2d', 'interp3d', 'lookup']

def lookup(x, y, z, table, mod='column'):
    if mod == 'column':
        table = table.T
    _x, _y, _z, _fs = table[0], table[1], table[2], table[3:]
    result = array([x, y, z])
    for _f in _fs:
        _table = array([_x, _y, _z, _f]).T
	result.append(trunc(interp3d(x, y, z, _table), 14))
    return result
"""    
filename='/home/vraphael/.python/spectre/lib/tsmc/crn40lp/layout_effect.dat'
from rawdata import *
layout_effect = base('layout_effect.dat').read()
t = layout_effect[0]
core = array([t['Lr'], t['Wr'], t['Nr'], t['SPBA'], t['SPA2'], t['SPA3'], t['SPA1'], t['SCA'], t['SCC'], t['SCB'], t['SPA'], t['SAP'], t['SODX1'], t['SODX2'], t['SPBA1'], t['SA3'], t['SA2'], t['SA1'], t['SODX'], t['SA6'], t['SA5'], t['SA4'], t['SAPB'], t['SODY'], t['SA']])
print lookup(0.041, 2, 16, core, mod='row')
"""
