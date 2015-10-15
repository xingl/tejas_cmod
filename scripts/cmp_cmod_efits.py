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
import math

e = 1.6*10**(-19)
a = 2.6863997038399379E-01

parser = op.OptionParser()
options,args = parser.parse_args()
lores_efit_file_name = args[0]
hires_efit_file_name = args[1]

### read EFIT file
psip_n_lores, Rgrid_lores, Zgrid_lores, F_lores, p_lores, ffprime_lores, pprime_lores, psirz_lores, qpsi_lores, rmag_lores, zmag_lores, nw_lores, psiax_lores, psisep_lores, Bctr_lores = read_EFIT_file(lores_efit_file_name)

psip_n_hires, Rgrid_hires, Zgrid_hires, F_hires, p_hires, ffprime_hires, pprime_hires, psirz_hires, qpsi_hires, rmag_hires, zmag_hires, nw_hires, psiax_hires, psisep_hires, Bctr_hires = read_EFIT_file(hires_efit_file_name)

rho_tor_spl_lores, rhot_n_lores = calc_rho_tor(psip_n_lores, psiax_lores, psisep_lores, qpsi_lores, nw_lores)

uni_rhot_lores = np.linspace(rhot_n_lores[0],rhot_n_lores[-1],1000)
q_unirhot_lores = interp(rhot_n_lores,qpsi_lores,uni_rhot_lores)
shat_lores = fd_d1_o4(q_unirhot_lores,uni_rhot_lores)*uni_rhot_lores/q_unirhot_lores

rho_tor_spl_hires, rhot_n_hires = calc_rho_tor(psip_n_hires, psiax_hires, psisep_hires, qpsi_hires, nw_hires)

uni_rhot_hires = np.linspace(rhot_n_hires[0],rhot_n_hires[-1],1000)
q_unirhot_hires = interp(rhot_n_hires,qpsi_hires,uni_rhot_hires)
shat_hires = fd_d1_o4(q_unirhot_hires,uni_rhot_hires)*uni_rhot_hires/q_unirhot_hires

plt.plot(uni_rhot_lores,shat_lores,'x',label='low res')
plt.plot(uni_rhot_hires,shat_hires,'r.',label='high res')
plt.xlabel('rho_tor_n')
plt.ylabel('shat')
plt.legend()
plt.show()



plt.plot(psip_n_lores,abs(qpsi_lores),'x',label='low res')
plt.plot(psip_n_hires,abs(qpsi_hires),'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('q')
plt.legend()
plt.show()

plt.plot(psip_n_lores,F_lores,'x',label='low res')
plt.plot(psip_n_hires,F_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('F')
plt.legend()
plt.show()


plt.plot(psip_n_lores,pprime_lores,'x',label='low res')
plt.plot(psip_n_hires,pprime_hires,'.',label='high res')
#plt.plot(psip_n_hires,2*math.pi*pprime_hires,'.',label='high resx2*pi')
plt.xlabel('psi_pol_n')
plt.ylabel('pprime')
plt.legend()
plt.show()

plt.plot(psip_n_lores,ffprime_lores,'x',label='low res')
plt.plot(psip_n_hires,ffprime_hires,'.',label='high res')
#plt.plot(psip_n_hires,2*math.pi*ffprime_hires,'.',label='high resx2*pi')
plt.xlabel('psi_pol_n')
plt.ylabel('ffprime')
plt.legend()
plt.show()

plt.plot(psip_n_lores,F_lores,'x',label='low res')
plt.plot(psip_n_hires,F_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('F')
plt.legend()
plt.show()

print 'R of magnetic axis', rmag_lores, rmag_hires, abs(abs(rmag_lores)-abs(rmag_hires))*2/(abs(rmag_lores)+abs(rmag_hires))
print 'Z of magnetic axis', zmag_lores, zmag_hires, abs(abs(zmag_lores)-abs(zmag_hires))*2/(abs(zmag_lores)+abs(zmag_hires))
print 'psi at magnetic axis', psiax_lores, psiax_hires
print 'psi at separatrix', psisep_lores, psisep_hires 
print 'difference of psisep-psiax', abs(abs(psisep_lores-psiax_lores)-abs(psisep_hires-psiax_hires))*2/(abs(psisep_lores-psiax_lores)+abs(psisep_hires-psiax_hires))
print 'Bctr', Bctr_lores, Bctr_hires, abs(Bctr_lores-Bctr_hires)*2/(Bctr_lores+Bctr_hires)

rdim_lores,zdim_lores,rctr_lores,rmin_lores,zmid_lores,current_lores = read_EFIT_parameters(lores_efit_file_name)

rdim_hires,zdim_hires,rctr_hires,rmin_hires,zmid_hires,current_hires = read_EFIT_parameters(hires_efit_file_name)

#print 'rdim', rdim_lores, rdim_hires, abs(rdim_lores-rdim_hires)*2/(rdim_lores+rdim_hires)
#print 'zdim', zdim_lores, zdim_hires, abs(zdim_lores-zdim_hires)*2/(zdim_lores+zdim_hires)
#print 'rctr', rctr_lores, rctr_hires, abs(rctr_lores-rctr_hires)*2/(rctr_lores+rctr_hires)
#print 'rmin', rmin_lores, rmin_hires, abs(rmin_lores-rmin_hires)*2/(rmin_lores+rmin_hires)
#print 'zmid', zmid_lores, zmid_hires, abs(zmid_lores-zmid_hires)*2/(zmid_lores+zmid_hires)
#print 'current', current_lores, current_hires, abs(current_lores-current_hires)*2/(current_lores+current_hires)



psip_n_obmp_lores, R_obmp_lores, B_pol_lores, B_tor_lores = calc_B_fields(Rgrid_lores, rmag_lores, Zgrid_lores, zmag_lores, psirz_lores, psiax_lores, psisep_lores, F_lores, nw_lores, psip_n_lores)

p_obmp_lores = interp(psip_n_lores,p_lores,psip_n_obmp_lores)
beta_obmp_lores = 8*math.pi*p_obmp_lores/(B_pol_lores**2+B_tor_lores**2)

psip_n_obmp_hires, R_obmp_hires, B_pol_hires, B_tor_hires = calc_B_fields(Rgrid_hires, rmag_hires, Zgrid_hires, zmag_hires, psirz_hires, psiax_hires, psisep_hires, F_hires, nw_hires, psip_n_hires)

p_obmp_hires = interp(psip_n_hires,p_hires,psip_n_obmp_hires)
beta_obmp_hires = 8*math.pi*p_obmp_hires/(B_pol_hires**2+B_tor_hires**2)

plt.plot(psip_n_obmp_lores,B_pol_lores,'x',label='low res')
plt.plot(psip_n_obmp_hires,B_pol_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('B_pol')
plt.legend()
plt.show()

plt.plot(psip_n_obmp_lores,B_tor_lores,'x',label='low res')
plt.plot(psip_n_obmp_hires,B_tor_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('B_tor')
plt.legend()
plt.show()


plt.plot(psip_n_obmp_lores,beta_obmp_lores,'x',label='low res')
plt.plot(psip_n_obmp_hires,beta_obmp_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('beta')
plt.legend()
plt.show()

pprime_obmp_lores = interp(psip_n_lores,pprime_lores,psip_n_obmp_lores)
ffprime_obmp_lores = interp(psip_n_lores,ffprime_lores,psip_n_obmp_lores)
jt_obmp_lores = R_obmp_lores*pprime_obmp_lores+ffprime_obmp_lores/R_obmp_lores

pprime_obmp_hires = interp(psip_n_hires,pprime_hires,psip_n_obmp_hires)
ffprime_obmp_hires = interp(psip_n_hires,ffprime_hires,psip_n_obmp_hires)
jt_obmp_hires = R_obmp_hires*pprime_obmp_hires+ffprime_obmp_hires/R_obmp_hires

#plt.plot(psip_n_obmp_lores,pprime_obmp_lores,'x',label='low res')
#plt.plot(psip_n_obmp_hires,pprime_obmp_hires,'r.',label='high res')
#plt.ylabel('pprime')
#plt.xlabel('psi_pol_n')
#plt.legend()
#plt.show()

#plt.plot(psip_n_obmp_lores,ffprime_obmp_lores,'x',label='low res')
#plt.plot(psip_n_obmp_hires,ffprime_obmp_hires,'r.',label='high res')
#plt.ylabel('ffprime')
#plt.xlabel('psi_pol_n')
#plt.legend()
#plt.show()

plt.plot(psip_n_obmp_lores,jt_obmp_lores,'x',label='low res')
plt.plot(psip_n_obmp_hires,jt_obmp_hires,'.',label='high res')
#plt.plot(psip_n_obmp_hires,2*math.pi*jt_obmp_hires,'r.',label='high resx2*pi')
plt.xlabel('psi_pol_n')
plt.ylabel('jt')
plt.legend()
plt.show()


#binfo_lores = np.genfromtxt('Binfo_g1120907032.01012')
#plt.plot(psip_n_obmp_lores,B_tor_lores,'x',label='x')
#plt.plot(binfo_lores[:,1],binfo_lores[:,3],'r.')
#plt.legend()
#plt.xlabel('psi_pol_n')
#plt.ylabel('B_tor')
#plt.show()

#plt.plot(psip_n_obmp_lores,B_pol_lores,'x',label='x')
#plt.plot(binfo_lores[:,1],binfo_lores[:,2],'r.')
#plt.legend()
#plt.xlabel('psi_pol_n')
#plt.ylabel('B_pol')
#plt.show()
#
#binfo_hires = np.genfromtxt('Binfo_g_901_901_1415_q')
#plt.plot(psip_n_obmp_hires,B_tor_hires,'x',label='x')
#plt.plot(binfo_hires[:,1],binfo_hires[:,3],'r.')
#plt.legend()
#plt.xlabel('psi_pol_n')
#plt.ylabel('B_tor')
#plt.show()

#plt.plot(psip_n_obmp_hires,B_pol_hires,'x',label='x')
#plt.plot(binfo_hires[:,1],binfo_hires[:,2],'r.')
#plt.legend()
#plt.xlabel('psi_pol_n')
#plt.ylabel('B_pol')
#plt.show()
