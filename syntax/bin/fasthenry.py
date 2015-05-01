from simulator import fasthenry as simulator
from syntax import *
from math import pi, sqrt
from simulator import fasthenry as simulator


__all__ = ['FastHenry', 'Title', 'Units', 'Default', 'Node', 'Segment', 'Port', 'Freq', 'Equiv']

class FastHenry(list):
    def __init__(self):
        list.__init__(self, [])
    def __str__(self):
        self.append('.end')
        return "\n".join([str(line) for line in self])
    def run(self, verbose=False, **options):
        simu = simulator(self)
        simu.run(verbose=verbose, **options)
        self.raw = simu.raw
    @property
    def stdout(self):
        s = []
        n = len(self.raw)
        nfreq = len(self.raw['freq'])
        for ifreq in xrange(nfreq):
            freq = self.raw['freq'][ifreq]
            if freq>=1e9:
                s.append( 'Frequency = %g GHz' % (freq/1e9) )
            elif freq>=1e6:
                s.append( 'Frequency = %g MHz' % (freq/1e6) )
            elif freq>=1e3:
                s.append( 'Frequency = %g kHz' % (freq/1e3) )
            else:
                s.append( 'Frequency = %g Hz' % (freq) )
            if n == 2:
                z11 = self.raw['z11'][ifreq]
                s.append( '  r11= %.4g Ohms'%(z11.real) )
                s.append( '  l11= %.4g nH'%(z11.imag/freq/2/pi*1e9) )
            elif n == 5:
                z11 = self.raw['z11'][ifreq]
                z12 = self.raw['z12'][ifreq]
                z22 = self.raw['z22'][ifreq]
                s.append( '  r11= %.4g Ohms'%(z11.real) )
                s.append( '  l11= %.4g nH'%(z11.imag/freq*1e9/2/pi) )
                s.append( '  m12= %.4g'%(z12.imag/sqrt(z11.imag*z22.imag) ) )
                s.append( '  r22= %.4g Ohms'%(z22.real) )
                s.append( '  l22= %.4g nH'%(z22.imag/freq*1e9/2/pi) )
                
            s.append('')
        return "\n".join(s)

