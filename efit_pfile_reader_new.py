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

e = 1.6*10**(-19)
a = 2.69E-01
impurity_charge = 5.0
ion_charge = 1.0

### three flags: -e (plot EFIT-file), -p (plot p-file), -s (shift Ti profile)
parser = op.OptionParser()
parser.add_option('--plot_EFIT','-e',action='store_const',const=1,help='Plot pressure,q,shat,Bt,Bp from EFIT file',default=0)
parser.add_option('--plot_pfile','-p',action='store_const',const=1,help='Plot ti,te,ni,ne from p-file',default=0)
parser.add_option('--shift_Ti','-s',action='store_const',const=1,help='shift Ti profile out by 0.005',default=0)
### take one argument: 'Imode' or 'Hmode'
options,args = parser.parse_args()
mode_type = args[0]
if mode_type == 'Imode' :
    efit_file_name = 'g1120907032.01012'
    p_file_name = 'p1120907032.01012'
elif mode_type == 'Hmode' :
    efit_file_name = 'g1120815027.01075_001'
    p_file_name = 'p1120815027.01075_001'
plot_EFIT = options.plot_EFIT
plot_pfile = options.plot_pfile
shift_Ti = options.shift_Ti

### read EFIT file
psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)

### find conversion relation from psip_n to rhot_n (normalized sqrt of psit)
### rho_tor_spl takes psi_pol (not normalized) and calculates rhot (normalized)
rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)

### output q vs rhot_n
#f = open('q_profile.dat','w')
#f.write('#1.rhot_n 2.q\n')
#np.savetxt(f,np.column_stack((rhot_n,qpsi)))
#f.close()

if plot_EFIT:
    ### plot pressure vs rhot_n
    plt.plot(rhot_n,p*1E-3,'x',label='EFIT')
    plt.xlabel('rhot_n')
    plt.ylabel('Pressure (kPa)')
    #plt.legend()
    plt.title(efit_file_name)
    plt.show()

    ### plot q vs rhot_n
    plt.plot(rhot_n,qpsi,'x',label='EFIT')
    plt.xlabel('rhot_n')
    plt.ylabel('q')
    #plt.legend()
    plt.title(efit_file_name)
    plt.show()

### calculate shat and plot it vs rhot_n
uni_rhot = np.linspace(rhot_n[0],rhot_n[-1],1000)
q_unirhot = interp(rhot_n,qpsi,uni_rhot)
shat = fd_d1_o4(q_unirhot,uni_rhot)*uni_rhot/q_unirhot

### output shat vs rhot_n
#f = open('shat_profile.dat','w')
#f.write('#1.rhot_n 2.shat\n')
#np.savetxt(f,np.column_stack((uni_rhot,shat)))
#f.close()

#if plot_EFIT:
# plot shat vs uni_rhot
#    plt.plot(uni_rhot,shat,'r.',label='EFIT')
#    plt.xlabel('rhot_n')
#    plt.ylabel('shat')
#    #plt.legend()
#    plt.title(efit_file_name)
#    plt.show()

### calculate B_pol and B_tor at the outboard midplane (obmp) uniform R grid
psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F, nw, psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

if plot_EFIT:
    ### plot B_pol vs rhot_n
    plt.plot(rhot_n_obmp,B_pol,'r.',label='EFIT')
    plt.xlabel('rhot_n')
    plt.ylabel('B_pol')
    #plt.legend()
    plt.title(efit_file_name)
    plt.show()

    ### plot B_tor vs rhot_n
    plt.plot(rhot_n_obmp,B_tor,'r.',label='EFIT')
    plt.xlabel('rhot_n')
    plt.ylabel('B_tor')
    #plt.legend()
    plt.title(efit_file_name)
    plt.show()

### read from p-file, raw data with each quantity having its own grid
### ne, ni are in the unit of 10^20 m^(-3)
### te, ti are in the unit of keV
psipne, ne, psipte, te, psipni, ni, psipti, ti = read_cmod_pfile_raw(p_file_name)
if plot_pfile:
    ### plot Te (keV) vs psip_te
    plt.plot(psipte,te,'x',label='Te p-file')
    ### plot Ti (keV) vs psip_ti
    plt.plot(psipti,ti,'x',label='Ti p-file')
    plt.xlabel('psip')
    plt.ylabel('keV')
    plt.legend()
    x1,x2,y1,y2=plt.axis()
    plt.axis((0.95,1.0,y1,y2))
    plt.title('raw '+p_file_name)
    plt.show()

    ### plot ne (10^20 m^(-3)) vs psip_ne
    plt.plot(psipne,ne,'x',label='ne p-file')
    ### plot ni (10^20 m^(-3)) vs psip_ni
    plt.plot(psipni,ni,'x',label='ni p-file')
    plt.xlabel('psip')
    plt.ylabel('10^20 m^-3')
    plt.legend()
    plt.title('raw '+p_file_name)
    x1,x2,y1,y2=plt.axis()
    plt.axis((0.95,1.0,y1,y2))
    plt.show()


### read pfile and calculate impurity density
if shift_Ti :
    psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name,shift_Ti=True,shift_Te=False,add_impurity=True)
else:
    psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name,shift_Ti=False,shift_Te=False,add_impurity=True)

psip = psi0*(psisep-psiax)+psiax
rhot0 = rho_tor_spl(psip)

### convert to SI units to calculate pressure later
ne = ne*1E20
te = te*e*1E03
ni = ni*1E20
ti = ti*e*1E03
nz = nz*1E20

### Zeff = (ni+nz*z^2)/(ni+nz*z)
zeff = (ni*ion_charge**2+nz*impurity_charge**2)/ne
### output Zeff vs rhot
#f = open('Zeff_profile.dat','w')
#f.write('#1.rhot_n 2.Zeff\n')
#np.savetxt(f,np.column_stack((rhot0,zeff)))
#f.close()
plt.plot(rhot0,zeff*te/ti)
plt.axvline(x=0.985)
plt.axvline(x=0.980)
plt.axvline(x=0.975)
plt.xlabel('rhot_n')
plt.ylabel('tau=Zeff*Te/Ti')
plt.show()
print "check quasineutrality:", max(abs(ne-(ni+impurity_charge*nz))/ne[0])

if plot_pfile:
    # plot Ti&Te(keV) vs rhot_n
    plt.plot(rhot0,ti/e*1E-03,'x',label='Ti')
    plt.plot(rhot0,te/e*1E-03,'x',label='Te')
    plt.xlabel('rhot_n')
    plt.ylabel('KeV')
    plt.legend()
    plt.axis((0.95,1.0,0,1))
    plt.title(p_file_name)
    plt.show()

    # plot ni,ne&nz(10^20 m^-3) vs rhot_n
    plt.plot(rhot0,ni*1E-20,'x',label='ni')
    plt.plot(rhot0,ne*1E-20,'x',label='ne')
    plt.plot(rhot0,nz*1E-20,'x',label='nz')
    plt.xlabel('rhot_n')
    plt.ylabel('10^20 m^-3')
    plt.legend()
    plt.axis((0.95,1.0,0,1.5))
    plt.title(p_file_name)
    plt.show()

    # plot Zeff vs rhot_n
    plt.plot(rhot0,zeff,'x')
    plt.xlabel('rhot_n')
    plt.ylabel('Zeff')
    plt.legend()
    plt.title(p_file_name)
    x1,x2,y1,y2=plt.axis()
    plt.axis((0.95,1.0,y1,y2))
    plt.show()

### plot p vs psip_n
plt.plot(rhot0,(ni*ti+ne*te+nz*ti)*1E-03,'x',label='p-file')
plt.plot(rhot_n,p*1E-03,'r.',label='EFIT file')
plt.xlabel('rhot_n')
plt.ylabel('Pressure (kPa)')
plt.axis((0.95,1.0,0,45))
#plt.title(mode_type)
plt.legend()
plt.show()

# output p-file info
#f = open('p_info.dat','w')
#f.write('# 1.rhot_n 2.ne 3.te 4.ni 5.ti 6.nz 7.psip_n\n # \n')
#np.savetxt(f,np.column_stack((rhot0,ne,te,ni,ti,nz,psi0)))
#f.close()
