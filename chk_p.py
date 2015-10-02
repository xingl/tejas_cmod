from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
import re
from interp import *
from finite_differences import *
from read_EFIT_file import *
from read_cmod_pfile import *
from calc_gammaE import *
from read_Er import *

e = 1.6*10**(-19)
impurity_charge = 5.0

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'

### read EFIT file
psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file_name)

#psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name,shift=True,add_impurity=True)
psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name,shift_Te=True,shift_Ti=True,add_impurity=True)
ne = ne*1E20
te = te*e*1E03
ni = ni*1E20
ti = ti*e*1E03
nz = nz*1E20

plt.plot(ti,te)
plt.show()

print "check quasineutrality:", max(abs(ne-(ni+impurity_charge*nz))/ne[0])
#plt.plot(psi0,(ni+impurity_charge*nz),'x',label='ni+z*nz')
#plt.plot(psi0,ne,'r.',label='ne')
#plt.legend()
#plt.show()

print "pressure different at boundary:", p[-1], (ni[-1]*ti[-1]+ne[-1]*te[-1]+nz[-1]*ti[-1])
dp = abs(p[-1]-(ni[-1]*ti[-1]+ne[-1]*te[-1]+nz[-1]*ti[-1]))
plt.plot(psi0,(ni*ti+ne*te+nz*ti),'x',label='p-file')
plt.plot(psip_n,p,'.',label='EFIT file')
plt.plot(psip_n,p+dp,'.',label='up-shifted EFIT file')
#plt.title('Te profile shifted to the right')
#plt.title('Te profile no shift')
plt.title('Ti&Te shifted')
plt.xlabel('rho_tor')
plt.ylabel('Pressure(Pa)')
plt.axis((0.9,1.0,0,40000))
plt.legend()
plt.show()
