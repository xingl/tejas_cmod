from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
from interp import *
from finite_differences import *
from calc_gammaE import *
from read_EFIT_file import *

parser = op.OptionParser()
options,args = parser.parse_args()
rbsProfs_file = args[0]
efit_file_name = args[1]

e = 1.6*10**(-19)
M = 3.3*10**(-27)

# read from rbsProfs
f=open(rbsProfs_file,'r')
rbs = f.read()
f.close()

rbs = rbs.split('\n')
a_factor = float(rbs[1].split()[3])
psip_sep_rbs = float(rbs[1].split()[2])
print "a_factor",a_factor
print "psisep",psip_sep_rbs

rbs = np.genfromtxt(rbsProfs_file)
isep = np.argmin(abs(rbs[:,0]-1.0))
a = rbs[isep,22]
print "minor radius",a_factor*a

minor_a = a_factor*a

rhot_n_rbs = rbs[:,0]
psip_n_rbs = rbs[:,1]
p_rbs = rbs[:,2]
n = rbs[:,3]*10**20
Ti = rbs[:,4]*e*10**3
Te = rbs[:,5]*e*10**3
gamE_N = rbs[:,9]
gamE0_N = rbs[:,10]
Er = rbs[:,16]
Er0 = rbs[:,17]
R_rbs = rbs[:,24]
omega_rbs = rbs[:,27]
Bp_rbs = rbs[:,25]
Bt_rbs = rbs[:,26]

unif_R = linspace(R_rbs[0],R_rbs[-1],len(R_rbs))
rhot_unifR = interp(R_rbs,rhot_n_rbs,unif_R)
p_unifR_rbs = interp(rhot_n_rbs,p_rbs,rhot_unifR)
n_unifR = interp(rhot_n_rbs,n,rhot_unifR)
Ti_unifR = interp(rhot_n_rbs,Ti,rhot_unifR)
Te_unifR = interp(rhot_n_rbs,Te,rhot_unifR)
Er_unifR = interp(rhot_n_rbs,Er,rhot_unifR)
Er0_unifR = interp(rhot_n_rbs,Er0,rhot_unifR)
Bp_unifR = interp(rhot_n_rbs,Bp_rbs,rhot_unifR)
Bt_unifR = interp(rhot_n_rbs,Bt_rbs,rhot_unifR)
psip_n_unifR = interp(R_rbs,psip_n_rbs,unif_R)

# read from EFIT file
psip_n,Rgrid,Zgrid,F,p,ffprime,pprime,psirz,qpsi,rmag,zmag,nw,psiax,psisep = read_EFIT_file(efit_file_name)

rho_tor_spl, rhot_n = calc_rho_tor(psip_n,psiax,psisep,qpsi,nw)

psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid,rmag,Zgrid,zmag,psirz,psiax,psisep,F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

B_pol_unifR = interp(R_obmp,B_pol,unif_R)
B_tor_unifR = interp(R_obmp,B_tor,unif_R)
p_unifR = interp(psip_n,p,psip_n_unifR)

gammaE_p = calc_gammaE_p(unif_R,n_unifR,Ti_unifR,p_unifR,B_tor_unifR,B_pol_unifR,minor_a)
gammaE_niti = calc_gammaE_niti(unif_R,n_unifR,Ti_unifR,n_unifR,Te_unifR,B_tor_unifR,B_pol_unifR,minor_a)
gammaE_Er = calc_gammaE_Er(unif_R,Ti_unifR,B_tor_unifR,B_pol_unifR,Er_unifR,minor_a)
gammaE_Er0 = calc_gammaE_Er(unif_R,Ti_unifR,B_tor_unifR,B_pol_unifR,Er0_unifR,minor_a)

Ti_obmp = interp(R_rbs,Ti,R_obmp)
Er_obmp = interp(R_rbs,Er,R_obmp)
plt.plot(R_obmp,Er_obmp,'x')
plt.plot(R_rbs,Er,'r.')
plt.show()
Er0_obmp = interp(R_rbs,Er0,R_obmp)
gammaE_hb,omega_t,rhot_unifpsip = calc_gammaE_hb(unif_R,psip_n_unifR*abs(psisep-psiax),Ti_unifR,B_tor_unifR,B_pol_unifR,Er_unifR,rhot_unifR,minor_a)
gammaE_hb0,omega_t0,rhot_unifpsip0 = calc_gammaE_hb(unif_R,psip_n_unifR*abs(psisep-psiax),Ti_unifR,B_tor_unifR,B_pol_unifR,Er0_unifR,rhot_unifR,minor_a)

plt.plot(rhot_unifpsip0,omega_t0,'x',label='hb0')
plt.plot(rhot_unifpsip,omega_t,'r.',label='hb')
plt.plot(rhot_n_rbs,omega_rbs,'s',label='rbs')
plt.legend()
plt.ylabel('omega_t')
plt.xlabel('rhot_n')
plt.show()

#plt.plot(rhot_unifR,abs(gammaE_niti),'x',label='niti')
#plt.plot(rhot_unifR,abs(gammaE_p),'.',label='p')
plt.plot(rhot_unifR,abs(gammaE_Er),'s',label='Er')
#plt.plot(rhot_unifR,abs(gammaE_Er0),'.',label='Er0')
plt.plot(rhot_unifpsip, abs(gammaE_hb),'^',label='hb')
#plt.plot(rhot_unifpsip0,abs(gammaE_hb0),'o',label='hb0')
plt.plot(rhot_n_rbs,abs(gamE_N),'.',label='rbs')
#plt.plot(rhot_n_rbs,abs(gamE0_N),'.',label='rbs0')
plt.ylabel('abs(gammaE)/(v_thi/a)')
plt.title(rbsProfs_file)
plt.legend()
#x1,x2,y1,y2 = plt.axis()
#plt.axis((0.95,1.0,0.0,100.0))
plt.show()

