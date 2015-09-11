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
psip_n_lores, Rgrid_lores, Zgrid_lores, F_lores, p_lores, ffprime_lores, pprime_lores, psirz_lores, qpsi_lores, rmag_lores, zmag_lores, nw_lores, psiax_lores, psisep_lores = read_EFIT_file(lores_efit_file_name)

psip_n_hires, Rgrid_hires, Zgrid_hires, F_hires, p_hires, ffprime_hires, pprime_hires, psirz_hires, qpsi_hires, rmag_hires, zmag_hires, nw_hires, psiax_hires, psisep_hires = read_EFIT_file(hires_efit_file_name)

plt.plot(psip_n_lores,qpsi_lores,'x',label='low res')
plt.plot(psip_n_hires/0.999,-qpsi_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('q')
plt.legend()
plt.show()

plt.plot(psip_n_lores,pprime_lores,'x',label='low res')
plt.plot(psip_n_hires,2*math.pi*pprime_hires,'r.',label='high res')

plt.xlabel('psi_pol_n')
plt.ylabel('pprime')
plt.legend()
plt.show()


plt.plot(psip_n_lores,ffprime_lores,'x',label='low res')
plt.plot(psip_n_hires,2*math.pi*ffprime_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('ffprime')
plt.legend()
plt.show()

psip_n_obmp_lores, R_obmp_lores, B_pol_lores, B_tor_lores = calc_B_fields(Rgrid_lores, rmag_lores, Zgrid_lores, zmag_lores, psirz_lores, psiax_lores, psisep_lores, F_lores, nw_lores, psip_n_lores)

p_obmp_lores = interp(psip_n_lores,p_lores,psip_n_obmp_lores)
beta_obmp_lores = 8*math.pi*p_obmp_lores/(B_pol_lores**2+B_tor_lores**2)

psip_n_obmp_hires, R_obmp_hires, B_pol_hires, B_tor_hires = calc_B_fields(Rgrid_hires, rmag_hires, Zgrid_hires, zmag_hires, psirz_hires, psiax_hires, psisep_hires, F_hires, nw_hires, psip_n_hires)

p_obmp_hires = interp(psip_n_hires,p_hires,psip_n_obmp_hires)
beta_obmp_hires = 8*math.pi*p_obmp_hires/(B_pol_hires**2+B_tor_hires**2)

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
plt.plot(psip_n_obmp_hires,2*math.pi*jt_obmp_hires,'r.',label='high res')
plt.xlabel('psi_pol_n')
plt.ylabel('jt')
plt.legend()
plt.show()



