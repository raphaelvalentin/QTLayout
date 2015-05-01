#from CSTlib import *


material = list()
material.append( Material(name='Aluminum', kappa=3.8528221e+07, tc1=0.00385, color=(0.75, 0.75, 0.75)) )
material.append( Material(name='Poly', kappa=0.0873362e+05, tc1=0.00385, color=(0.25, 1.0, 0.50)) )
material.append( Material(name='CopperTM2', kappa=5.9321833e+07, tc1=0.00385) )
material.append( Material(name='CopperM2M8', kappa=5.1020408e+07, tc1=0.00385) )
material.append( Material(name='CopperM1', kappa=4.6264168e+07, tc1=0.00385) )
material.append( Material(name='Oxide', epsilon=4.2, transparency=0.5, color=(1, 0.5, 0) ) )
material.append( Material(name='Vacuum', transparency=0.5, color=(0, 1, 1)) )
material.append( Material(name='Si', epsilon=11.9, kappa=10.0, transparency=0.5) )
material.append( Material(name='SiN', epsilon=7.0, transparency=0.5) )
material.append( Material(name='TEOS', epsilon=4.0, transparency=0.5) )
material.append( Material(name='USG', epsilon=4.0, transparency=0.5) )
material.append( Material(name='BD', epsilon=3.0, transparency=0.5) )
material.append( Material(name='SiC', epsilon=5.1, transparency=0.5) )
material.append( Material(name='SiO2', epsilon=4.0, transparency=0.5) )
material.append( Material(name='SiON', epsilon=5.0, transparency=0.5) )


Layer.unit = 'A'
dielectric = Layers()
dielectric.append( Layer(name='STI', material='SiO2', thickness=3100, zmax=0, eps=4) )
dielectric.append( Layer(name='ILD1a', material='SiN', thickness=350, zmin=0, eps=7) )     
dielectric.append( Layer(name='ILD1b', material='SiO2', thickness=2800, zmin=350, eps=4) )     
dielectric.append( Layer(name='IMD1a', material='SiC', thickness=300, zmin=3150, eps=5.1) )     
dielectric.append( Layer(name='IMD1b', material='BD', thickness=1150, zmin=3450, eps=3) )     
dielectric.append( Layer(name='IMD2a', material='SiC', thickness=500, zmin=4600, eps=5.1) )     
dielectric.append( Layer(name='IMD2b', material='BD', thickness=3050, zmin=5100, eps=3) )     
dielectric.append( Layer(name='IMD3a', material='SiC', thickness=500, zmin=8150, eps=5.1) )     
dielectric.append( Layer(name='IMD3b', material='BD', thickness=3050, zmin=8650, eps=3) )     
dielectric.append( Layer(name='IMD4a', material='SiC', thickness=500, zmin=11700, eps=5.1) )     
dielectric.append( Layer(name='IMD4b', material='BD', thickness=3050, zmin=12200, eps=3) )     
dielectric.append( Layer(name='IMD5a', material='SiC', thickness=500, zmin=15250, eps=5.1) )     
dielectric.append( Layer(name='IMD5b', material='BD', thickness=3050, zmin=15750, eps=3) )     
dielectric.append( Layer(name='IMD6a', material='SiC', thickness=500, zmin=18800, eps=5.1) )     
dielectric.append( Layer(name='IMD6b', material='BD', thickness=3050, zmin=19300, eps=3) )     
dielectric.append( Layer(name='UTVa', material='SiN', thickness=750, zmin=22350, eps=7) )     
dielectric.append( Layer(name='UTVb', material='USG', thickness=7000, zmin=23100, eps=4) )     
dielectric.append( Layer(name='UTMa', material='SiN', thickness=1500, zmin=30100, eps=7) )     
dielectric.append( Layer(name='UTMb', material='USG', thickness=32000, zmin=31600, eps=4) )     
dielectric.append( Layer(name='Pass1', material='SiN', thickness=750, zmin=63600, eps=7) )      
dielectric.append( Layer(name='Pass2', material='TEOS', thickness=4000, zmin=64350, eps=4) )     
dielectric.append( Layer(name='Pass3', material='SiN', thickness=750, zmin=68350, eps=7) )
dielectric.append( Layer(name='Pass4', material='TEOS', thickness=2500, zmin=69100, eps=4) )
dielectric.append( Layer(name='Pass5', material='TEOS', thickness=4000, zmin=71600, eps=4) )
dielectric.append( Layer(name='Pass6', material='SiN', thickness=6000, zmin=75600, eps=7) )


# calculate equivalent (electrical) thickness
eps = 3.0
Layer.unit = 'um'
dielectric['STI']   = Layer(name='STI',    zmax=0, thickness=eps/dielectric['STI'].eps*dielectric['STI'].thickness)
dielectric['ILD1a'] = Layer(name='ILD1a',  zmin=0, thickness=eps/dielectric['ILD1a'].eps*dielectric['ILD1a'].thickness)
dielectric['ILD1b'] = Layer(name='ILD1b',  zmin=dielectric['ILD1a'].zmax,  thickness=eps/dielectric['ILD1b'].eps*dielectric['ILD1b'].thickness)
dielectric['IMD1a'] = Layer(name='IMD1a',  zmin=dielectric['ILD1b'].zmax,  thickness=eps/dielectric['IMD1a'].eps*dielectric['IMD1a'].thickness)
dielectric['IMD1b'] = Layer(name='IMD1b',  zmin=dielectric['IMD1a'].zmax,  thickness=eps/dielectric['IMD1b'].eps*dielectric['IMD1b'].thickness)
dielectric['IMD2a'] = Layer(name='IMD2a',  zmin=dielectric['IMD1b'].zmax,  thickness=eps/dielectric['IMD2a'].eps*dielectric['IMD2a'].thickness)
dielectric['IMD2b'] = Layer(name='IMD2b',  zmin=dielectric['IMD2a'].zmax,  thickness=eps/dielectric['IMD2b'].eps*dielectric['IMD2b'].thickness)
dielectric['IMD3a'] = Layer(name='IMD3a',  zmin=dielectric['IMD2b'].zmax,  thickness=eps/dielectric['IMD3a'].eps*dielectric['IMD3a'].thickness)
dielectric['IMD3b'] = Layer(name='IMD3b',  zmin=dielectric['IMD3a'].zmax,  thickness=eps/dielectric['IMD3b'].eps*dielectric['IMD3b'].thickness)
dielectric['IMD4a'] = Layer(name='IMD4a',  zmin=dielectric['IMD3b'].zmax,  thickness=eps/dielectric['IMD4a'].eps*dielectric['IMD4a'].thickness)
dielectric['IMD4b'] = Layer(name='IMD4b',  zmin=dielectric['IMD4a'].zmax,  thickness=eps/dielectric['IMD4b'].eps*dielectric['IMD4b'].thickness)
dielectric['IMD5a'] = Layer(name='IMD5a',  zmin=dielectric['IMD4b'].zmax,  thickness=eps/dielectric['IMD5a'].eps*dielectric['IMD5a'].thickness)
dielectric['IMD5b'] = Layer(name='IMD5b',  zmin=dielectric['IMD5a'].zmax,  thickness=eps/dielectric['IMD5b'].eps*dielectric['IMD5b'].thickness)
dielectric['IMD6a'] = Layer(name='IMD6a',  zmin=dielectric['IMD5b'].zmax,  thickness=eps/dielectric['IMD6a'].eps*dielectric['IMD6a'].thickness)
dielectric['IMD6b'] = Layer(name='IMD6b',  zmin=dielectric['IMD6a'].zmax,  thickness=eps/dielectric['IMD6b'].eps*dielectric['IMD6b'].thickness)
  
eps = 4.0
dielectric['UTVa'] = Layer(name='UTVa',  zmin=dielectric['IMD6b'].zmax,  thickness=eps/dielectric['UTVa'].eps*dielectric['UTVa'].thickness)
dielectric['UTVb'] = Layer(name='UTVb',  zmin=dielectric['UTVa'].zmax,  thickness=eps/dielectric['UTVb'].eps*dielectric['UTVb'].thickness)
dielectric['UTMa'] = Layer(name='UTMa',  zmin=dielectric['UTVb'].zmax,  thickness=eps/dielectric['UTMa'].eps*dielectric['UTMa'].thickness)
dielectric['UTMb'] = Layer(name='UTMb',  zmin=dielectric['UTMa'].zmax,  thickness=eps/dielectric['UTMb'].eps*dielectric['UTMb'].thickness)
dielectric['Pass1'] = Layer(name='Pass1',  zmin=dielectric['UTMb'].zmax,  thickness=eps/dielectric['Pass1'].eps*dielectric['Pass1'].thickness)
dielectric['Pass2'] = Layer(name='Pass2',  zmin=dielectric['Pass1'].zmax,  thickness=eps/dielectric['Pass2'].eps*dielectric['Pass2'].thickness)
dielectric['Pass3'] = Layer(name='Pass3',  zmin=dielectric['Pass2'].zmax,  thickness=eps/dielectric['Pass3'].eps*dielectric['Pass3'].thickness)
dielectric['Pass4'] = Layer(name='Pass4',  zmin=dielectric['Pass3'].zmax,  thickness=eps/dielectric['Pass4'].eps*dielectric['Pass4'].thickness)
dielectric['Pass5'] = Layer(name='Pass5',  zmin=dielectric['Pass4'].zmax,  thickness=eps/dielectric['Pass5'].eps*dielectric['Pass5'].thickness)
dielectric['Pass6'] = Layer(name='Pass6',  zmin=dielectric['Pass5'].zmax,  thickness=eps/dielectric['Pass6'].eps*dielectric['Pass6'].thickness)
    
Layer.unit = 'um'
metal = Layers()
metal.append( Layer(name='ALRDL', material='Aluminum', zmin=dielectric['Pass5'].zmin, thickness=1.45) )
metal.append( Layer(name='TM2', material='CopperTM2', zmin=dielectric['UTMa'].zmin, thickness=3.4) )
metal.append( Layer(name='Metal6', material='CopperM2M8', zmax=dielectric['IMD6b'].zmax, thickness=0.2) )
metal.append( Layer(name='Metal5', material='CopperM2M8', zmax=dielectric['IMD5b'].zmax, thickness=0.2) )
metal.append( Layer(name='Metal4', material='CopperM2M8', zmax=dielectric['IMD4b'].zmax, thickness=0.2) )
metal.append( Layer(name='Metal3', material='CopperM2M8', zmax=dielectric['IMD3b'].zmax, thickness=0.2) )
metal.append( Layer(name='Metal2', material='CopperM2M8', zmax=dielectric['IMD2b'].zmax, thickness=0.2) )
metal.append( Layer(name='Metal1', material='CopperM1', zmax=dielectric['IMD1b'].zmax, thickness=0.165) )
metal.append( Layer(name='Poly', material='Poly', zmin=0.035, thickness=0.1) )
metal.append( Layer(name='AA', material='Si', zmax=0.0, thickness=0.31) )
    
material.append( Material(name='Oxide1', epsilon=3, transparency=0.5) )
material.append( Material(name='Oxide2', epsilon=4, transparency=0.5) )

# equivalent IMD
z0 = -750.0
z1 = dielectric['STI'].zmin
z2 = dielectric['IMD6b'].zmax
z3 = metal['ALRDL'].zmax+0.5
z4 = 750.0
Layer.unit = 'um'
dielectric = Layers()
dielectric.append( Layer(name='PSUB', material='Si', zmin=z0, zmax=z1) )
dielectric.append( Layer(name='IMD1', material='Oxide1', zmin=z1, zmax=z2) )
dielectric.append( Layer(name='IMD2', material='Oxide2', zmin=z2, zmax=z3) )
dielectric.append( Layer(name='Air', material='Vacuum', zmin=z3, zmax=z4) )
    


        

