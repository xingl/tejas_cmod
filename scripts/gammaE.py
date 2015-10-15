#! /usr/bin/python

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

parser = op.OptionParser(description='calculate gamma_ExB')
parser.add_option('--plot_gammaE','-p',action='store_const',const=1,help='Plot gammaE vs rhot_n',default=1)
parser.add_option('--shift_Ti','-s',action='store_const',const=1,help='shift Ti profile out by 0.005',default=1)
parser.add_option('--plot_lingr','-l',action='store_const',const=1,help='Plot gammaE/linear_growth_rate vs rhot_n',default=0)
parser.add_option('--test_calc','-t',action='store_const',const=1,help='test calculaitons of gammaE',default=0)
options,args = parser.parse_args()
mode_type = args[0]

if mode_type == 'Imode' :
    efit_file_name = 'g_901_901_1415_q'
    p_file_name = 'p1120907032.01012'
    Er_file_name = 'Er1120907032.01010_v20140623'
elif mode_type == 'Hmode' :
    efit_file_name = 'g1120815027.01075_001'
    p_file_name = 'p1120815027.01075_001'
else:
    exit("""
Please include EFIT, p-file, Er-file(if available) as arguments (in that order).
     """)
plot_gammaE = options.plot_gammaE
shift_Ti = options.shift_Ti
plot_lingr = options.plot_lingr
test_calc = options.test_calc

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax,psisep = read_EFIT_file(efit_file_name)

rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)

psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

p_obmp = interp(psip_n,p,psip_n_obmp)

if shift_Ti :
    psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name)
else:
    psi0, ne, te, ni, ti, nz=read_cmod_pfile(p_file_name,shift_Ti=False,shift_Te=False)

ne = ne*1E20
te = te*e*1E03
ni = ni*1E20
ti = ti*e*1E03
nz = nz*1E20

rhot0=rho_tor_spl(psi0*(psisep-psiax)+psiax)

ni_obmp = interp(psi0,ni,psip_n_obmp)
ti_obmp = interp(psi0,ti,psip_n_obmp)
ne_obmp = interp(psi0,ne,psip_n_obmp)
te_obmp = interp(psi0,te,psip_n_obmp)

if mode_type == 'Imode' :
    if shift_Ti:
        shift_Er = True
    else:
        shift_Er = False
    Er_shift = 0.005  #Shift in poloidal flux coord
    psi0_Er, Er, Er_error = read_Er(Er_file_name,shift_Er,Er_shift)
    rhot0_Er = interp(psip_n_obmp,rhot_n_obmp,psi0_Er)
    if test_calc:
        plt.plot(rhot0_Er, Er)
        plt.xlabel('rhot_n')
        plt.ylabel('Er(kV/m)')
        plt.title(mode_type)
        plt.show()

    Er = Er*1E03
    # psi0_Er starts from 0.9 to 1.055, psip_obmp goes from 0 to 1.06 
    # high end of psip_n_obmp has to be outside of high end of psi0_Er
    # conversion makes Er outside the original range not useful  
    psi_Er_f = np.argmin(abs(psip_n_obmp-psi0_Er[0]))
    if psi0_Er[-1] <= psip_n_obmp[-1] :
    	psi_Er_l = np.argmin(abs(psip_n_obmp-psi0_Er[-1]))
    else:
        psi_Er_l = np.argmin(abs(psip_n_obmp-1.01))
    Er_obmp = interp(psi0_Er,Er,psip_n_obmp[psi_Er_f:psi_Er_l+1])
    gammaE_Er = calc_gammaE_Er(R_obmp[psi_Er_f:psi_Er_l+1],te_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,a)
    gammaE_hb, omega_t, rhot_unipsip = calc_gammaE_hb(R_obmp[psi_Er_f:psi_Er_l+1],psip_obmp[psi_Er_f:psi_Er_l+1],te_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp,rhot_n_obmp[psi_Er_f:psi_Er_l+1],a)
elif mode_type == 'Hmode' :
    psi_Er_f = np.argmin(abs(psip_n_obmp-0.947))
    psi_Er_l = np.argmin(abs(psip_n_obmp-1.01))

Er_gradp,gammaE_niti = calc_gammaE_niti(R_obmp[psi_Er_f:psi_Er_l+1],ni_obmp[psi_Er_f:psi_Er_l+1],ti_obmp[psi_Er_f:psi_Er_l+1],te_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1],B_pol[psi_Er_f:psi_Er_l+1],a)
gammaE_hb_gradp, omega_t, rhot_unipsip2 = calc_gammaE_hb(R_obmp[psi_Er_f:psi_Er_l+1],psip_obmp[psi_Er_f:psi_Er_l+1],te_obmp[psi_Er_f:psi_Er_l+1],B_tor[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_gradp,rhot_n_obmp[psi_Er_f:psi_Er_l+1],a)

if test_calc:
    plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_niti),'.',label='gammaE_niti')
    plt.plot(rhot_unipsip2, abs(gammaE_hb_gradp),'^', label='hb_gradp')
    if mode_type == 'Imode' :
        plt.plot(rhot_unipsip, abs(gammaE_hb),'^', label='hb')
        plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_Er),'x',label='exp Er')
    plt.xlabel('rho_tor')
    plt.ylabel('gammaE/(c_s/a)')
    plt.legend()
    plt.title(mode_type)
    plt.axis((0.95,1.0,0,55))
    plt.show()

if plot_gammaE and plot_lingr:
    plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_niti)/2.5,'.',label='gammaE_niti/2.5')
    if mode_type == 'Imode' :
        plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],abs(gammaE_Er)/2.5,'x',label='gammaE_Er/2.5')
    plt.xlabel('rhot_n')
    plt.ylabel('c_s/a')
    plt.legend()
    plt.title(mode_type)
    #plt.axis((0.95,1.0,0,55))
    plt.axis((0.95,1.0,0,25))
    #plt.show()

    #plt.plot(rhot_n_obmp[psi_Er_f:psi_Er_l+1],ni_obmp[psi_Er_f:psi_Er_l+1]*ti_obmp[psi_Er_f:psi_Er_l+1],'.',label='ni*ti')
    #plt.xlabel('rhot_n')
    #plt.ylabel('P_i(Pa)')
    #plt.title(mode_type)
    #plt.axis((0.95,1.0,0,20000))
    #plt.legend()
    #plt.show()
    x0_values = [92,96,965,97,975,98,985,99,995]
    gamE_gamL = np.empty(len(x0_values))
    gamNT_gamL = np.empty(len(x0_values))
    x=np.empty(len(x0_values))
    x_gammaE = rhot_n_obmp[psi_Er_f:psi_Er_l+1]
    x_gammaE_fine = linspace(x_gammaE[0],x_gammaE[-1],20*len(x_gammaE))
    gamE = interp(x_gammaE,gammaE_Er,x_gammaE_fine)
    gamNT = interp(x_gammaE,gammaE_niti,x_gammaE_fine)
    for i in range(len(x0_values)):
        filename = 'scan'+str(x0_values[i])+'.dat'
        scandat = np.genfromtxt(filename)
        gr = scandat[:,4]
        x0 = scandat[:,1]
        plt.scatter(x0,gr,alpha=0.5,color='red')
        x[i] = np.mean(x0[:])
        x0_ind = np.argmin(abs(x_gammaE_fine-x0[0]))
        gamE_gamL[i] = float(abs(gamE[x0_ind])/max(gr))
        gamNT_gamL[i] = float(abs(gamNT[x0_ind])/max(gr))
    plt.show()

    plt.scatter(x,gamNT_gamL,marker='s',color='blue',label='gammaE_niti/max(gammaL)')
    plt.scatter(x,gamE_gamL,marker='o',color='red',label='gammaE_Er/max(gammaL)')
    plt.xlabel('rhot_n')
    plt.ylabel('ratio')
    plt.title(mode_type)
    plt.axis((0.95,1.0,0,20))
    plt.legend()
    plt.show()
    
    f = open('gamE_gamL.dat','w')
    f.write('#1.rhot_n 2.gamE/gamL 3.gamNT/gamL \n#\n')
    np.savetxt(f,np.column_stack((x,gamE_gamL,gamNT_gamL)))
    f.close()
