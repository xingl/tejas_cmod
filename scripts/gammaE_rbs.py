#! /usr/bin/python

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
parser.add_option('--test_calc_gammaE','-t',action='store_const',const=1,help='test calc_gammaE subroutines using Bt,Bp from rbsProfs',default=0)
options,args = parser.parse_args()
rbsProfs_file = args[0]
efit_file_name = args[1]
test_calc_gammaE = options.test_calc_gammaE

e = 1.6*10**(-19)
M = 3.3*10**(-27)

f=open(rbsProfs_file,'r')
rbs = f.read()
f.close()

rbs = rbs.split('\n')
a_factor = float(rbs[1].split()[3])
psip_sep_rbs = float(rbs[1].split()[2])

rbs = np.genfromtxt(rbsProfs_file)
isep = np.argmin(abs(rbs[:,0]-1.0))
a = rbs[isep,22]

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

if test_calc_gammaE:
    Er_gradp_rbs,gammaE_rbs_niti = calc_gammaE_niti(unif_R,n_unifR,Ti_unifR,Te_unifR,Bt_unifR,Bp_unifR,minor_a)
    gammaE_rbs_Er = calc_gammaE_Er(unif_R,Te_unifR,Bt_unifR,Bp_unifR,Er_unifR,minor_a)
    gammaE_rbs_Er0 = calc_gammaE_Er(unif_R,Te_unifR,Bt_unifR,Bp_unifR,Er0_unifR,minor_a)
    gammaE_rbs_hb,omega_t,rhot_unifpsip = calc_gammaE_hb(R_rbs,psip_n_rbs*abs(psip_sep_rbs),Te,Bt_rbs,Bp_rbs,Er,rhot_n_rbs,minor_a)
    gammaE_rbs_hb0,omega_t0,rhot_unifpsip0 = calc_gammaE_hb(R_rbs,psip_n_rbs*abs(psip_sep_rbs),Te,Bt_rbs,Bp_rbs,Er0,rhot_n_rbs,minor_a)

    plt.plot(rhot_unifR,abs(gammaE_rbs_niti),'^',label='niti')
    plt.plot(rhot_unifR,abs(gammaE_rbs_Er0),'x',label='Er0')
    plt.plot(rhot_unifpsip0,abs(gammaE_rbs_hb0),'.',label='hb0')
    plt.plot(rhot_n_rbs,abs(gamE0_N),'.',label='rbs0')
    plt.ylabel('abs(gammaE0_rbs)/(c_s/a)')
    plt.title(rbsProfs_file)
    plt.legend()
    plt.axis((0.95,1.0,0.0,50.0))
    plt.show()

    plt.plot(rhot_unifR,abs(gammaE_rbs_Er),'s',label='Er')
    plt.plot(rhot_unifpsip,abs(gammaE_rbs_hb),'x',label='hb')
    plt.plot(rhot_n_rbs,abs(gamE_N),'.',label='rbs')
    plt.ylabel('abs(gammaE_rbs)/(c_s/a)')
    plt.title(rbsProfs_file)
    plt.legend()
    plt.axis((0.95,1.0,0.0,50.0))
    plt.show()

psip_n,Rgrid,Zgrid,F,p,ffprime,pprime,psirz,qpsi,rmag,zmag,nw,psiax,psisep = read_EFIT_file(efit_file_name)

rho_tor_spl, rhot_n = calc_rho_tor(psip_n,psiax,psisep,qpsi,nw)

psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid,rmag,Zgrid,zmag,psirz,psiax,psisep,F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

B_pol_unifR = interp(R_obmp,B_pol,unif_R)
B_tor_unifR = interp(R_obmp,B_tor,unif_R)
p_unifR = interp(psip_n,p,psip_n_unifR)

Er_gradp,gammaE_niti = calc_gammaE_niti(unif_R,n_unifR,Ti_unifR,Te_unifR,B_tor_unifR,B_pol_unifR,minor_a)

gammaE_Er = calc_gammaE_Er(unif_R,Te_unifR,B_tor_unifR,B_pol_unifR,Er_unifR,minor_a)

Ti_obmp = interp(R_rbs,Ti,R_obmp)
Te_obmp = interp(R_rbs,Te,R_obmp)
Er_obmp = interp(R_rbs,Er,R_obmp)

ped_f = np.argmin(abs(psip_n_obmp-0.947))
ped_l = np.argmin(abs(psip_n_obmp-0.997))
#gammaE_hb,omega_t,rhot_unifpsip = calc_gammaE_hb(R_obmp[ped_f:ped_l+1],psip_n_obmp[ped_f:ped_l+1]*abs(psisep),Te_obmp[ped_f:ped_l+1],B_tor[ped_f:ped_l+1],B_pol[ped_f:ped_l+1],Er_obmp[ped_f:ped_l+1],rhot_n_obmp[ped_f:ped_l+1],minor_a)

plt.plot(rhot_unifR,abs(gammaE_niti),'x',label='niti')
plt.plot(rhot_unifR,abs(gammaE_Er),'.',label='Er')
plt.xlabel('rhot_n')
plt.ylabel('abs(gammaE)/(c_s/a)')
plt.title(rbsProfs_file)
plt.legend()
x1,x2,y1,y2 = plt.axis()
plt.axis((0.95,1.0,0.0,100.0))
plt.show()

