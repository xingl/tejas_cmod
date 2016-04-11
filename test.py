from read_EFIT_file import *
from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *

#efit_file_name = 'g_new_901_901_1415'
#efit_file_name = 'FS2RZ_gEFIT'
#flux_surface = '0.97'
parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = args[0]
flux_surface = args[1]

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)
print 'psiax = ', psiax
print 'psisep = ', psisep
fs=float(flux_surface)

rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)

psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F, nw, psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

print len(Rgrid), len(Zgrid)
print rmag,zmag

ngrid_r = 3
ngrid_z = 3

if (1==1):
    Z_fs = np.empty(0) 
    R_fs_out_grid = np.empty(0) 
    B_pol_fs_out = np.empty(0)
    B_tor_fs_out = np.empty(0)
    R_fs_in_grid = np.empty(0)
    B_pol_fs_in = np.empty(0)
    B_tor_fs_in = np.empty(0)

    Z0_ind = np.argmin(abs(Zgrid-zmag))
    print Z0_ind

    #decrease Z from midplane until the lowest psip at that Z is larger than fs 
    while (np.min((psirz[Z0_ind,:]-psiax)/(psisep-psiax))<fs):

        psi_pol = psirz[Z0_ind,:]
	#plt.plot(Rgrid,psi_pol)
	#plt.show()
        #Rmag_ind is the position of lowest psip at that Z
        #Rmag_ind = np.argmin(abs(psi_pol))
        #print 'Z0_ind = ', Z0_ind
	#print 'Z = ', Zgrid[Z0_ind]
        ###Rmag_ind = np.argmin(psi_pol)
        Rmag_ind = np.argmin(abs(psi_pol))
	#print 'Rmag_ind = ', Rmag_ind
	#print 'Rmag = ', Rgrid[Rmag_ind]
	#R_fs_ind_out is the index of psip_n_temp that is on fs
        psi_pol_out = psi_pol[Rmag_ind:].copy()
        psip_n_temp = (psi_pol_out-psiax)/(psisep-psiax)
        R_fs_ind_out = np.argmin(abs(psip_n_temp-fs))
	#unif_R_out is a fine grid of Rgrid on the outer
        unif_R_out = np.linspace(Rgrid[Rmag_ind],Rgrid[-1],nw*10)
        psip_n_unifR_out = interp(Rgrid[Rmag_ind:],psip_n_temp,unif_R_out)
        #psifs_ind_out is the postition of psip_n = fs
        psifs_ind_out = np.argmin(abs(psip_n_unifR_out-fs))
        #psip_n_fs_out is the local grid of psip_n around fs 
        psip_n_fs_out = psip_n_unifR_out[psifs_ind_out-ngrid_r:psifs_ind_out+ngrid_r].copy()
        #R_fs_out is the local grid of R around fs 
        R_fs_out = unif_R_out[psifs_ind_out-ngrid_r:psifs_ind_out+ngrid_r].copy()
        #B_pol_Z_out is the z component of B_pol at fs, B_pol_Z = 1/R d psip /d R
        B_pol_Z_out = fd_d1_o4(psip_n_fs_out*(psisep-psiax)+psiax,R_fs_out)/R_fs_out
	#B_tor = F/R
        F_out = interp(psip_n,F,psip_n_fs_out)
        B_tor_out = F_out/R_fs_out
     
        #psi_pol_z selects ngrid_z points above and below (R_fs,Z0) from psirz
        psi_pol_z = psirz[Z0_ind-ngrid_z:Z0_ind+ngrid_z,R_fs_ind_out+Rmag_ind]
	#unif_z is a local fine grid of Z
        unif_z = np.linspace(Zgrid[Z0_ind-ngrid_z],Zgrid[Z0_ind+ngrid_z],10*ngrid_z)
	#psi_pol_unifz is a local grid of psi_pol
        psi_pol_unifz = interp(Zgrid[Z0_ind-ngrid_z:Z0_ind+ngrid_z],psi_pol_z,unif_z)
        #B_pol_R_out is the R component of B_pol at fs, B_pol_R = 1/R d psip/d Z
        B_pol_R_out = fd_d1_o4(psi_pol_unifz,unif_z)/R_fs_out[ngrid_r]
        #z_fs_ind is the position of Z0_ind in the newly constructed array unif_z
        z_fs_ind=np.argmin(abs(unif_z-Zgrid[Z0_ind]))
        #B_pol_out is the total B_pol field at psip=fs
        B_pol_out = np.sqrt(B_pol_Z_out[ngrid_r]**2+B_pol_R_out[z_fs_ind]**2) 

        #similar procedure at inner side
        psi_pol_in = psi_pol[:Rmag_ind].copy()
        psip_n_temp = (psi_pol_in-psiax)/(psisep-psiax)
        R_fs_ind_in = np.argmin(abs(psip_n_temp-fs))
        unif_R_in = np.linspace(Rgrid[0],Rgrid[Rmag_ind],nw*10)
        psip_n_unifR_in = interp(Rgrid[:Rmag_ind],psip_n_temp,unif_R_in)
        psifs_ind_in = np.argmin(abs(psip_n_unifR_in-fs))
        psip_n_fs_in = psip_n_unifR_in[psifs_ind_in-ngrid_r:psifs_ind_in+ngrid_r].copy()
        R_fs_in = unif_R_in[psifs_ind_in-ngrid_r:psifs_ind_in+ngrid_r].copy()
        B_pol_Z_in = fd_d1_o4(psip_n_fs_in*(psisep-psiax)+psiax,R_fs_in)/R_fs_in
        F_in = interp(psip_n,F,psip_n_fs_in)
        B_tor_in = F_in/R_fs_in

        psi_pol_z = psirz[Z0_ind-ngrid_z:Z0_ind+ngrid_z,R_fs_ind_in]
        psi_pol_unifz = interp(Zgrid[Z0_ind-ngrid_z:Z0_ind+ngrid_z],psi_pol_z,unif_z)
        B_pol_R_in = fd_d1_o4(psi_pol_unifz,unif_z)/R_fs_in[ngrid_r]
        B_pol_in = np.sqrt(B_pol_Z_in[ngrid_r]**2+B_pol_R_in[z_fs_ind]**2) 

        Z_fs = np.append(Z_fs,Zgrid[Z0_ind])
        R_fs_out_grid = np.append(R_fs_out_grid,R_fs_out[ngrid_r])
        B_pol_fs_out = np.append(B_pol_fs_out,B_pol_out)
        B_tor_fs_out = np.append(B_tor_fs_out,B_tor_out[ngrid_r])
        R_fs_in_grid = np.append(R_fs_in_grid,R_fs_in[ngrid_r])
        B_pol_fs_in = np.append(B_pol_fs_in,B_pol_in)
        B_tor_fs_in = np.append(B_tor_fs_in,B_tor_in[ngrid_r])

        Z0_ind = Z0_ind-1

    Z_fs_r = np.flipud(Z_fs)
    R_fs_out_grid_r = np.flipud(R_fs_out_grid)
    B_pol_fs_out_r = np.flipud(B_pol_fs_out)
    B_tor_fs_out_r = np.flipud(B_tor_fs_out)

    B_tot_fs_out_r = np.sqrt(B_tor_fs_out_r**2+B_pol_fs_out_r**2)
    B_tot_fs_in = np.sqrt(B_tor_fs_in**2+B_pol_fs_in**2)

    theta_out_r = np.arctan(abs(Z_fs_r)/R_fs_out_grid_r)
    theta_in = -np.arctan(abs(Z_fs)/R_fs_in_grid)+np.pi

    R_fs = np.concatenate([R_fs_in_grid,R_fs_out_grid_r])
    B_pol = np.concatenate([B_pol_fs_in,B_pol_fs_out_r])
    B_tor = np.concatenate([B_tor_fs_in,B_tor_fs_out_r])
    B_tot = np.concatenate([B_tot_fs_in,B_tot_fs_out_r])
    theta_fs = np.concatenate([theta_in,theta_out_r])

plt.plot(R_fs,B_pol,label='B_pol')
plt.plot(R_fs,abs(B_tor),label='B_tor')
plt.plot(R_fs,B_tot,label='B_tot')
plt.title(efit_file_name+' at flux surface of '+flux_surface+' +- 0.01')
plt.xlabel('R')
plt.ylabel('B')
plt.legend()
plt.show()

f = open('Bp.dat','w')
f.write('# R Bp Bt\n # \n')
np.savetxt(f,np.column_stack((R_fs,B_pol,B_tor)))
f.close()
