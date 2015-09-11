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
a = 2.6863997038399379E-01

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'

### read EFIT file
psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax,psisep = read_EFIT_file(efit_file_name)

### find conversion relation from psi_pol to rhot_n
rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)
plt.plot(rhot_n,p)
plt.ylabel('p')
plt.xlabel('rho_tor')
#x1,x2,y1,y2=plt.axis()
#plt.axis((0.945,1.0,y1,y2))
plt.show()

### calculate B_pol and B_tor at the outboard midplane (obmp)
psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name)
ne = ne*1E20
te = te*e*1E03
ni = ni*1E20
ti = ti*e*1E03
nz = nz*1E20

#plt.plot(psi0,te,'x',label='electron temp')
#plt.plot(psi0,ti,'r.',label='ion temp')
#plt.legend()
#plt.show()

rhot0=rho_tor_spl(psi0*(psisep-psiax)+psiax)

ni_obmp = interp(psi0,ni,psip_n_obmp)
ti_obmp = interp(psi0,ti,psip_n_obmp)
ne_obmp = interp(psi0,ne,psip_n_obmp)
te_obmp = interp(psi0,te,psip_n_obmp)

Er_file = 'Er1120907032.01010_v20140623'
shift_Er = True   #Shift Er profile in poloidal flux coord
Er_shift = 0.005  #Shift in poloidal flux coord
psi0_Er, Er, Er_error = read_Er(Er_file,shift_Er,Er_shift)
Er = Er*1E03

# psi0_Er starts from 0.9 to 1.055, psip_obmp goes from 0 to 1.06 
# high end of psip_n_obmp has to be outside of high end of psi0_Er
# conversion makes Er outside the original range not useful  
psi_Er_f = np.argmin(abs(psip_n_obmp-psi0_Er[0]))
psi_Er_l = np.argmin(abs(psip_n_obmp-psi0_Er[-1]))
Er_obmp = interp(psi0_Er,Er,psip_n_obmp[psi_Er_f:psi_Er_l+1])

#plt.plot(psi0_Er,Er,'x')
#plt.plot(psip_n_obmp[psi_Er_f:psi_Er_l+1],Er_obmp,'r.')
#plt.ylabel('Er')
#plt.show()

p_obmp = interp(psip_n,p,psip_n_obmp)

gammaE_p = calc_gammaE_p(R_obmp,ni_obmp,ti_obmp,p_obmp,B_tor,B_pol,a)
gammaE_niti = calc_gammaE_niti(R_obmp,ni_obmp,ti_obmp,ne_obmp,te_obmp,B_tor,B_pol,a)
gammaE_Er = calc_gammaE_Er(R_obmp[psi_Er_f:psi_Er_l+1], ti_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,a)
gammaE_hb, omega_t, rhot_unipsip = calc_gammaE_hb(R_obmp[psi_Er_f:psi_Er_l+1],psip_obmp[psi_Er_f:psi_Er_l+1],ti_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,rhot_n_obmp[psi_Er_f:psi_Er_l+1],a)

#plt.plot(rhot_n_obmp,abs(gammaE_p),'x',label='EFIT P')
#plt.plot(rhot_n_obmp,abs(gammaE_niti),'.',label='exp ni*ti')
#plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_Er),'x',label='exp Er')
#plt.plot(rhot_unipsip, abs(gammaE_hb),'.', label='hb')
#plt.xlabel('rho_tor')
#plt.ylabel('gammaE/(v_thi/a)')
#plt.legend()
#x1,x2,y1,y2 = plt.axis()
#plt.axis((0.95,1.0,0.0,200.0))
#plt.show()


gammaE_hb_cs, omega_t, rhot_unipsip = calc_gammaE_hb(R_obmp[psi_Er_f:psi_Er_l+1],psip_obmp[psi_Er_f:psi_Er_l+1],te_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,rhot_n_obmp[psi_Er_f:psi_Er_l+1],a)

plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_Er)/2.5,'rx',label='gammaE_Er/2.5')
#x1,x2,y1,y2 = plt.axis()
#plt.axis((0.945,1.0,0.0,y2))
plt.xlabel('rho_tor')
plt.ylabel('gamma/(C_s/a)')
plt.legend()
#plt.show()

gr095 = [0.7,1.1,1.2,1.3,1.2,1.15,1.05]
x095 = np.repeat(0.95,len(gr095))
plt.scatter(x095,gr095,alpha=0.5)

gr096 = [0.55,0.9,1.0,1.1,1.0,0.95,0.85]
x096 = np.repeat(0.96,len(gr096))
plt.scatter(x096,gr096,alpha=0.5)

gr0965 = [0.4,0.7,0.8,0.85,0.8,0.75,0.6]
x0965 = np.repeat(0.965,len(gr0965))
plt.scatter(x0965,gr0965,alpha=0.5)

gr097 = [0.55,0.8,0.95,0.9,0.85,0.75,0.7]
x097 = np.repeat(0.97,len(gr097))
plt.scatter(x097,gr097,alpha=0.5)

gr098 = [1.0,1.1,1.5,2.0,2.5,3.0,3.5,4.0]
x098 = np.repeat(0.98,len(gr098))
plt.scatter(x098,gr098,alpha=0.5)

#data0985 = np.genfromtxt('scan.log')
#gr0985 = data0985[:,4]
gr0985 = [2,2.5,3,4,5,6]
x0985 = np.repeat(0.985,len(gr0985))
plt.scatter(x0985,gr0985,alpha=0.5)
x1,x2,y1,y2 = plt.axis()
plt.axis((0.945,1.0,0.0,y2))
plt.show()

