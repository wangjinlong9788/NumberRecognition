__author__ = 'wangjinlong'
# argument are inner radius, outer radius, epsilon, Nmodes, zmax
from pydefault import *
from DWFA_cyl_func_Ng import *
import os
#print sys.argv

# structure parameters:
# becareful Ng switch convention wrt Gai (a>b in Ng!!!)
b      = 1e-3#float(sys.argv[1])
a      = 1.3e-3#float(sys.argv[2])
epsilon= 4.41#float(sys.argv[3])
Nmode  = 1#int(float(sys.argv[4]))
zmax   = 1e-3#float(sys.argv[5])
#fin    = sys.argv[6]
#fout   = sys.argv[7]

#command='cp '+fin+' '+fout
#os.system(command)

print 'argument:'
print 'inner_radius =', b
print 'outer_radius =', a
print 'expsilon_rel =', epsilon
print 'NumberofModes=', Nmode
print 'Maximum z pos=', zmax

### initial set that works
#b = 450e-6
#a = 550e-6
#epsilon = 4.41 #dielectric constant of the medium
mu=1.0000
# number of mode we want to get
#Nmode = 1


# if fin sdds file is given zmax is OVERWRTTTEN ** need a condition here eventually **


cms=299792458.0
zmax= 0.00098790855413#cms*np.abs((max_t-min_t))

print 'Maximum z pos=', zmax

# compute the Green's function
Nz=10000
zmin=0#-zmax#0
# zmax is changed to have conservative headroom
zmax=2.*zmax

r=b/2.0
r0=b/2.0#1.e-3/2.0
azimuthal_mode=1  #m=0

RootAmplit, RootWavVec= FindMode(a,b,azimuthal_mode,mu,epsilon,Nmode,1)
n=azimuthal_mode
zGreen, WlGreen = Long_GreenFunction (RootAmplit, RootWavVec, r0, r, b, a, n, zmin, zmax, Nz, mu, epsilon)
zGreen, WtGreen = Trans_GreenFunction(RootAmplit, RootWavVec, r0, r, b, a, n, zmin, zmax, Nz, mu, epsilon)
plt.plot(zGreen, WtGreen)
plt.show()
