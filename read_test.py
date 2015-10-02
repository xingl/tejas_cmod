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
from read_iterdb_file import *
from w_iterdb import *

e = 1.6*10**(-19)
a = 2.6863997038399379E-01

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = 'g1120907032.01012'
p_file_name = 'p1120907032.01012'
iterdb_filename = args[0]

### read EFIT file
psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)

### find conversion relation from psi_pol to rhot_n
rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)

uni_rhot = np.linspace(rhot_n[0],rhot_n[-1],1000)
q_unirhot = interp(rhot_n,qpsi,uni_rhot)
shat = fd_d1_o4(q_unirhot,uni_rhot)*uni_rhot/q_unirhot

#plt.plot(rhot_n,qpsi,'x')
#plt.plot(uni_rhot,q_unirhot,'r.')
#plt.show()

#plt.plot(uni_rhot,shat)
#plt.xlabel('rho_tor_n')
#plt.ylabel('shat')
#plt.show()

### calculate B_pol and B_tor at the outboard midplane (obmp) uniform R grid
psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

#plt.plot(psip_n_obmp,rhot_n_obmp,'x')
#plt.plot(psip_n,rhot_n,'r.')
#plt.show()

### read from p-file
### ne, ni are in the unit of 10^20 m^(-3)
### te, ti are in the unit of KeV
psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name)
ne = ne*1E20
te = te*e*1E03
ni = ni*1E20
ti = ti*e*1E03
nz = nz*1E20

##plt.plot(psi0,te,'x',label='electron temp')
##plt.plot(psi0,ti,'r.',label='ion temp')
##plt.legend()
##plt.show()

rhot0=rho_tor_spl(psi0*(psisep-psiax)+psiax)

ni_obmp = interp(psi0,ni,psip_n_obmp)
ti_obmp = interp(psi0,ti,psip_n_obmp)
ne_obmp = interp(psi0,ne,psip_n_obmp)
te_obmp = interp(psi0,te,psip_n_obmp)

rhotp_obmp = interp(psi0,rhot0,psip_n_obmp)
q_obmp = interp(psip_n, qpsi, psip_n_obmp)

#plt.plot(psip_n_obmp,rhot_n_obmp,'x')
#plt.plot(psi0,rhot0,'r.')
#plt.show()

### read from Er file
### Er is in the unit of kV/m
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

##plt.plot(psi0_Er,Er,'x')
##plt.plot(psip_n_obmp[psi_Er_f:psi_Er_l+1],Er_obmp,'r.')
##plt.ylabel('Er')
##plt.show()

p_obmp = interp(psip_n,p,psip_n_obmp)

#gammaE_p = calc_gammaE_p(R_obmp,ni_obmp,ti_obmp,p_obmp,B_tor,B_pol,a)
#gammaE_niti = calc_gammaE_niti(R_obmp,ni_obmp,ti_obmp,ne_obmp,te_obmp,B_tor,B_pol,a)
#gammaE_Er = calc_gammaE_Er(R_obmp[psi_Er_f:psi_Er_l+1], ti_obmp[psi_Er_f:psi_Er_l+1],ti_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,a)
#gammaE_hb, omega_t, rhot_unipsip = calc_gammaE_hb(R_obmp[psi_Er_f:psi_Er_l+1],psip_obmp[psi_Er_f:psi_Er_l+1],ti_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,rhot_n_obmp[psi_Er_f:psi_Er_l+1],a)

##plt.plot(rhot_n_obmp,abs(gammaE_p),'x',label='EFIT P')
##plt.plot(rhot_n_obmp,abs(gammaE_niti),'.',label='exp ni*ti')
#plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_Er),'x',label='exp Er')
#plt.plot(rhot_unipsip, abs(gammaE_hb),'.', label='hb')
#plt.xlabel('rho_tor_n')
#plt.ylabel('gammaE/(v_thi/a)')
#plt.legend()
#plt.show()


###plt.plot(psip_n_obmp[psi_Er_f:psi_Er_l+1],rhot_n_obmp[psi_Er_f:psi_Er_l+1],'x')
###plt.plot(psi0,rhot0,'r.')
###plt.show()

#gammaE_hb_cs, omega_t, rhot_unipsip = calc_gammaE_hb(R_obmp[psi_Er_f:psi_Er_l+1],psip_obmp[psi_Er_f:psi_Er_l+1],te_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,rhot_n_obmp[psi_Er_f:psi_Er_l+1],a)

##plt.plot(rhot_unipsip, abs(gammaE_hb_cs)/2.5, '.', label ='gammaE_hb_cs/2.5')
##x1,x2,y1,y2 = plt.axis()
##plt.axis((0.95,1.0,y1,y2))
##plt.legend()
##plt.show()

### read from iterdb file
### te, ti are in the unit of ev
### ne, ni are in the unit of m^(-3)
rhot_idb, te_idb, ti_idb, ne_idb, ni_idb, nb_idb, vrot_idb = read_iterdb_file(iterdb_filename)

psip_idb = interp(rhot_n,psip_n,rhot_idb)
psipne, ne, psipte, te, psipni, ni, psipti, ti = read_cmod_pfile_raw(p_file_name)

plt.plot(psipte,te*1E03,'x',label='te p-file')
plt.plot(psip_idb,te_idb,'r.',label='te iterdb')
plt.xlabel('rhot_n')
plt.legend()
plt.show()

plt.plot(psipti,ti*1E03,'x',label='ti p-file')
plt.plot(psip_idb,ti_idb,'r.',label='ti iterdb')
plt.xlabel('rhot_n')
plt.legend()
plt.show()

plt.plot(psipne,ne*1E20,'x',label='ne p-file')
plt.plot(psip_idb,ne_idb,'r.',label='ne iterdb')
plt.xlabel('rhot_n')
plt.legend()
plt.show()

plt.plot(psipni,ni*1E20,'x',label='ni p-file')
plt.plot(psip_idb,ni_idb,'r.',label='ni iterdb')
plt.xlabel('rhot_n')
plt.legend()
plt.show()

#plt.plot(psip_n,rhot_n,'x')
#plt.plot(psip_n_obmp, rhot_n_obmp,'r.')
#plt.plot(psip_idb, rhot_idb,'s')
#plt.show()

#rhot_gene,gammaE_gene,omega_tor_gene = calc_gammaE_gene(R_obmp[psi_Er_f:psi_Er_l+1], rhot_n_obmp[psi_Er_f:psi_Er_l+1], te_obmp[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp, q_obmp[psi_Er_f:psi_Er_l+1], a)

#psip_gene = interp(rhot_n,psip_n,rhot_gene)
#plt.plot(psip_gene, omega_tor_gene,'x')
#plt.plot(psip_idb,vrot_idb,'r.')
#plt.ylabel('gammaE_gene')
#plt.xlabel('rhot_n')
#plt.show()


