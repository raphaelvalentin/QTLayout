from math import *

width = 25.
spacing = 15.


ldc=2.45e-11*15/width*(1-2.14e-11/2.45e-11*(1-width/15)) 
csi=2.025e-13*(1+0.02524*width+0.01698*spacing)/(1+0.02524*15+0.01698*10) 
cox=1.27e-13*(1+0.269*width+0.06075*spacing)/(1+0.269*15+0.06075*10) 
tc1=0.0039
ls=3.028e-10*(log(4900/(width+0.07291*spacing))-0.5)/(log(4900/(15+0.07291*10))-0.5) 
rdc=0.1928/width*15
cs=8e-15/spacing*10
fac=1e+10
coupling=0.4808*(log(563.2/(spacing+0.05922*width))-1)/(log(563.2/(10+0.05922*15))-1) 
rsi=64.78*(1+-0.008018*width+-0.003427*spacing)/(1+-0.008018*15+-0.003427*10) 
rac=0.5268/sqrt(width/15) 

print 'W%dS%d'%(width,spacing)
print '%.2e'%rdc
print '%.2e'%ls
print '%.2e'%rac
print '%.2e'%ldc
print '%.2e'%cs
print '%.2e'%cox
print '%.2e'%rsi
print '%.2e'%csi
print '%.2e'%coupling
