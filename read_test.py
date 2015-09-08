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

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = args[0]
p_file_name = args[1]

### read EFIT file
psi_pol, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmid, nw = read_EFIT_file(efit_file_name)

### find conversion relation from psi_pol to rho_tor 
rho_tor_spl, psi_pol, rho_tor = calc_rho_tor(psi_pol, qpsi, nw)
#plt.plot(psi_pol,rho_tor,'x')
#plt.xlabel('psi_pol')
#plt.ylabel('rho_tor')
#plt.show()

### calculate B_pol and B_tor at the outboard midplane (obmp)
psi_pol_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmid, psirz, psi_pol, F)
#plt.plot(psi_pol_obmp,'x',label='psi_pol at midplane')
#plt.plot(psi_pol,'r.',label='psi_pol grid of EFIT file fields')
#plt.legend()
#plt.show()

Binfo = np.genfromtxt('Binfo_g1120907032.01012')
plt.plot((psi_pol_obmp-psi_pol[0])/(psi_pol[-1]-psi_pol[0]), B_pol,'x')
plt.plot(Binfo[:,1],Binfo[:,2],'r.')
plt.title('B_pol')
plt.show()

plt.plot((psi_pol_obmp-psi_pol[0])/(psi_pol[-1]-psi_pol[0]), B_tor,'x')
plt.plot(Binfo[:,1],Binfo[:,3],'r.')
plt.title('B_tol')
plt.show()
### read pfile
psi0, ne, te, ni, ti=read_cmod_pfile(p_file_name)

rhot0=rho_tor_spl(psi0*(psi_pol[-1]-psi_pol[0])+psi_pol[0])
rhot0=rhot0/rhot0[-1]

eprof = np.genfromtxt('profiles_wb_gene_e.dat')
iprof = np.genfromtxt('profiles_wb_gene_i.dat')

plt.plot(eprof[:,0],eprof[:,2],'x')
plt.plot(rhot0, te,'r.')
plt.title('te')
plt.show()
plt.plot(eprof[:,0],eprof[:,3]/10.,'x')
plt.plot(rhot0, ne,'r.')
plt.title('ne')
plt.show()
plt.plot(iprof[:,0],iprof[:,2],'x')
plt.plot(rhot0, ti,'r.')
plt.title('ti')
plt.show()
plt.plot(iprof[:,0],iprof[:,3]/10.,'x')
plt.plot(rhot0, ni,'r.')
plt.title('ni')
plt.show()
