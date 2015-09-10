from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
from interp import *
from finite_differences import *
from calc_gammaE import *

parser = op.OptionParser()
options,args = parser.parse_args()
rbsProfs_file = args[0]

e = 1.6*10**(-19)
M = 3.3*10**(-27)

f=open(rbsProfs_file,'r')
rbs = f.read()
f.close()

rbs = rbs.split('\n')
a_factor = float(rbs[1].split()[3])
psip_sep = float(rbs[1].split()[2])
#print "a_factor",a_factor
#print "psip_sep",psip_sep

rbs = np.genfromtxt(rbsProfs_file)
isep = np.argmin(abs(rbs[:,0]-1.0))
#print "rhotor[isep]",rbs[isep,0]
a = rbs[isep,22]
#print "a[isep]",a
#print "a_min_eff",a_factor*a

minor_a = a_factor*a

rhot_n = rbs[:,0]
#plt.plot(rhot_n)
#plt.plot(rhot_n**2)
#plt.ylabel('rhot_n')
#plt.show()
psip_n = rbs[:,1]
#plt.plot(psip_n)
#plt.ylabel('psip_n')
#plt.show()
p = rbs[:,2]
n = rbs[:,3]*10**20
Ti = rbs[:,4]*e*10**3
Te = rbs[:,5]*e*10**3
gamE_N = rbs[:,9]
gamE0_N = rbs[:,10]
Er = rbs[:,16]
Er0 = rbs[:,17]
R = rbs[:,24]
omega = rbs[:,27]
#plt.plot(R)
#plt.ylabel('R')
#plt.show()
Bp = rbs[:,25]
Bt = rbs[:,26]

#plt.plot(rhot_n,p,'x',label='p')
#plt.plot(rhot_n,n*(Ti+Te),'r.',label='n*(Ti+Te)')
#plt.ylabel('Pressure(Pa)')
#plt.legend()
#plt.show()
unif_R = linspace(R[0],R[-1],len(R))
#plt.plot(R,'x')
#plt.plot(unif_R,'r.')
#plt.show()
#print "first and last R:", R[0], R[-1]
#print "first and last unif_R:", unif_R[0], unif_R[-1]
rhot_unifR = interp(R,rhot_n,unif_R)
#plt.plot(R,rhot_n,'x')
#plt.plot(unif_R,rhot_unifR,'r.')
#plt.show()
p_unifR = interp(rhot_n,p,rhot_unifR)
n_unifR = interp(rhot_n,n,rhot_unifR)
Ti_unifR = interp(rhot_n,Ti,rhot_unifR)
Te_unifR = interp(rhot_n,Te,rhot_unifR)
Er_unifR = interp(rhot_n,Er,rhot_unifR)
Er0_unifR = interp(rhot_n,Er0,rhot_unifR)
Bp_unifR = interp(rhot_n,Bp,rhot_unifR)
Bt_unifR = interp(rhot_n,Bt,rhot_unifR)


gammaE_p = calc_gammaE_p(unif_R,n_unifR,Ti_unifR,p_unifR,Bt_unifR,Bp_unifR,minor_a)

gammaE_niti = calc_gammaE_niti(unif_R,n_unifR,Ti_unifR,n_unifR,Te_unifR,Bt_unifR,Bp_unifR,minor_a)

gammaE_Er = calc_gammaE_Er(unif_R,Ti_unifR,Bt_unifR,Bp_unifR,Er_unifR,minor_a)
gammaE_Er0 = calc_gammaE_Er(unif_R,Ti_unifR,Bt_unifR,Bp_unifR,Er0_unifR,minor_a)

#psip_sep = np.trapz(unif_R*Bp_unifR,unif_R)
gammaE_hb,omega_t,rhot_unifpsip = calc_gammaE_hb(R,psip_n*abs(psip_sep),Ti,Bt,Bp,Er,rhot_n,minor_a)
gammaE_hb0,omega_t0,rhot_unifpsip0 = calc_gammaE_hb(R,psip_n*abs(psip_sep),Ti,Bt,Bp,Er0,rhot_n,minor_a)

#plt.plot(rhot_n,omega,'x',label='rbs')
#plt.plot(rhot_unifpsip,omega_t,'o',label='hb')
#plt.plot(rhot_unifpsip0,omega_t0,'^',label='hb0')
#plt.ylabel('Er/RBp')
#plt.legend()
#plt.show()

plt.plot(rhot_unifR,abs(gammaE_niti),'x',label='niti')
#plt.plot(rhot_n,abs(gamE0_N),'r.')
#plt.plot(rhot_unifR,abs(gammaE_p),'^',label='p')
plt.plot(rhot_unifR,abs(gammaE_Er),'s',label='Er')
plt.plot(rhot_unifR,abs(gammaE_Er0),'.',label='Er0')
plt.plot(rhot_unifpsip,abs(gammaE_hb),'.',label='hb')
plt.plot(rhot_unifpsip0,abs(gammaE_hb0),'.',label='hb0')
plt.plot(rhot_n,abs(gamE_N),'.',label='rbs')
plt.plot(rhot_n,abs(gamE0_N),'.',label='rbs0')
plt.ylabel('abs(gammaE)/(v_thi/a)')
plt.title(rbsProfs_file)
plt.legend()
plt.show()


