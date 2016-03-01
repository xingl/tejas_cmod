from pylab import *
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *
from read_EFIT_file import *

def fields_fs_upper(efit_file_name,flux_surface,ngrid_r,ngrid_z,extra_r,extra_z):
    
    fs = float(flux_surface)

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)


    Z_fs = np.empty(0) 
    R_fs_out_grid = np.empty(0) 
    B_pol_fs_out = np.empty(0)
    B_tor_fs_out = np.empty(0)
    R_fs_in_grid = np.empty(0)
    B_pol_fs_in = np.empty(0)
    B_tor_fs_in = np.empty(0)

    Z0_ind = np.argmin(abs(Zgrid-zmag))

    #decrease Z from midplane until the lowest psip at that Z is larger than fs 
    while (np.min((psirz[Z0_ind,:]-psiax)/(psisep-psiax))<fs):

        psi_pol = psirz[Z0_ind,:]
        #Rmag_ind is the position of lowest psip at that Z
        Rmag_ind = np.argmin(abs(psi_pol))
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


    extra_Rgrid = np.linspace(R_fs_in[ngrid_r],R_fs_out[ngrid_r],extra_r)
    extra_Zgrid = np.linspace(Zgrid[Z0_ind],Zgrid[Z0_ind+3],extra_z)
    Zgrid_frame = Zgrid[Z0_ind-ngrid_z:Z0_ind+1+ngrid_z]
    psip_extra_frame = np.empty([1+2*ngrid_z,extra_r])
    psip_extra = np.empty([extra_z,extra_r])
    for i in arange(0,2*ngrid_z+1):
	psi_pol = psirz[Z0_ind-ngrid_z+i,:]
    	psip_extra_frame[i,:] = interp(Rgrid,psi_pol,extra_Rgrid)
	#plt.plot(extra_Rgrid,psip_extra[i,:])
	#plt.show()
    for i in arange(0,extra_r):
	psip_z_frame = psip_extra_frame[:,i]
	psip_extra[:,i]=interp(Zgrid_frame,psip_z_frame,extra_Zgrid)
	#plt.plot(extra_Zgrid,psip_extra[:,i],'.')
	#plt.show()
    B_pol_bottom_out = np.empty(0)
    B_tor_bottom_out = np.empty(0)
    R_bottom_out = np.empty(0)
    B_pol_bottom_in = np.empty(0)
    B_tor_bottom_in = np.empty(0)
    R_bottom_in = np.empty(0)
    for i in arange(0,extra_z):
	psip_n_extra = (psip_extra[i,:]-psiax)/(psisep-psiax)
	r_ind = np.argmin(psip_n_extra)
	psip_n_extra_out = psip_n_extra[r_ind:].copy()
	extra_r_ind_out = np.argmin(psip_n_extra_out-fs)
	B_pol_Z_extra_out = \
	fd_d1_o4(psip_n_extra_out*(psisep-psiax)+psiax,extra_Rgrid[r_ind:])/extra_Rgrid[r_ind:]
	F_extra_out = interp(psip_n,F,psip_n_extra_out)
	B_tor_extra_out = F_extra_out/extra_Rgrid[r_ind:]
	B_pol_R_extra_out = \
	fd_d1_o4(psip_extra[:,r_ind+extra_r_ind_out],extra_Zgrid)/extra_Rgrid[r_ind+extra_r_ind_out]
	this_B_pol = np.sqrt(B_pol_R_extra_out[i]**2+B_pol_Z_extra_out[extra_r_ind_out]**2)
	if this_B_pol>1.E-3 :
		B_pol_bottom_out = np.append(B_pol_bottom_out,this_B_pol)
		R_bottom_out = np.append(R_bottom_out,extra_Rgrid[r_ind+extra_r_ind_out])
		B_tor_bottom_out = np.append(B_tor_bottom_out,B_tor_extra_out[extra_r_ind_out])
	
	psip_n_extra_in = psip_n_extra[:r_ind].copy()
	extra_r_ind_in = np.argmin(psip_n_extra_in-fs)
	B_pol_Z_extra_in = \
	fd_d1_o4(psip_n_extra_in*(psisep-psiax)+psiax,extra_Rgrid[:r_ind])/extra_Rgrid[:r_ind]
	F_extra_in = interp(psip_n,F,psip_n_extra_in)
	B_tor_extra_in = F_extra_in/extra_Rgrid[:r_ind]
	B_pol_R_extra_in = \
	fd_d1_o4(psip_extra[:,extra_r_ind_in],extra_Zgrid)/extra_Rgrid[extra_r_ind_in]
	this_B_pol = np.sqrt(B_pol_R_extra_in[i]**2+B_pol_Z_extra_in[extra_r_ind_in]**2)
	if this_B_pol>1.E-3 :
		B_pol_bottom_in = np.append(B_pol_bottom_in,this_B_pol)
		R_bottom_in = np.append(R_bottom_in,extra_Rgrid[extra_r_ind_in])
		B_tor_bottom_in = np.append(B_tor_bottom_in,B_tor_extra_in[extra_r_ind_in])

    R_bottom_out_r = np.flipud(R_bottom_out)
    B_pol_bottom_out_r = np.flipud(B_pol_bottom_out)
    B_tor_bottom_out_r = np.flipud(B_tor_bottom_out)
    B_tot_bottom_out_r = np.sqrt(B_tor_bottom_out_r**2+B_pol_bottom_out_r**2)
    B_tot_bottom_in = np.sqrt(B_tor_bottom_in**2+B_pol_bottom_in**2)

    R_fs_out_grid_r = np.flipud(R_fs_out_grid)
    B_pol_fs_out_r = np.flipud(B_pol_fs_out)
    B_tor_fs_out_r = np.flipud(B_tor_fs_out)
    B_tot_fs_out_r = np.sqrt(B_tor_fs_out_r**2+B_pol_fs_out_r**2)
    B_tot_fs_in = np.sqrt(B_tor_fs_in**2+B_pol_fs_in**2)

    #theta_out_r = np.arctan(abs(Z_fs)/R_fs_out_grid_r)
    #theta_in = -np.arctan(abs(Z_fs)/R_fs_in_grid)+np.pi

    R_fs = np.concatenate([R_fs_in_grid,R_bottom_in,R_bottom_out_r,R_fs_out_grid_r])
    B_pol = np.concatenate([B_pol_fs_in,B_pol_bottom_in,B_pol_bottom_out_r,B_pol_fs_out_r])
    B_tor = np.concatenate([B_tor_fs_in,B_tor_bottom_in,B_tor_bottom_out_r,B_tor_fs_out_r])
    B_tot = np.concatenate([B_tot_fs_in,B_tot_bottom_in,B_tot_bottom_out_r,B_tot_fs_out_r])
    #theta_fs = np.concatenate([theta_in,theta_out_r])

    #f=open('fs_'+efit_file_name+'.dat','w')
    #f.write('# R B_pol B_tor B_tot theta\n')
    #np.savetxt(f,np.column_stack((R_fs,B_pol,B_tor,B_tot,theta_fs)))
    #f.close()

    return R_fs, B_pol, B_tor, B_tot

def fields_fs_lower(efit_file_name,flux_surface,ngrid_r,ngrid_z):
    
    fs = float(flux_surface)

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)


    Z_fs = np.empty(0) 
    R_fs_out_grid = np.empty(0) 
    B_pol_fs_out = np.empty(0)
    B_tor_fs_out = np.empty(0)
    R_fs_in_grid = np.empty(0)
    B_pol_fs_in = np.empty(0)
    B_tor_fs_in = np.empty(0)

    Z0_ind = np.argmin(abs(Zgrid-zmag))

    #decrease Z from midplane until the lowest psip at that Z is larger than fs 
    while (np.min((psirz[Z0_ind,:]-psiax)/(psisep-psiax))<fs):

        psi_pol = psirz[Z0_ind,:]
        #Rmag_ind is the position of lowest psip at that Z
        Rmag_ind = np.argmin(abs(psi_pol))
        psi_pol_out = psi_pol[Rmag_ind:].copy()
        psip_n_temp = (psi_pol_out-psiax)/(psisep-psiax)
        R_fs_ind_out = np.argmin(abs(psip_n_temp-fs))
        #print ((psirz[Z0_ind,Rmag_ind+R_fs_ind_out]-psiax)/(psisep-psiax)-fs)/fs
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
        F_out = interp(psip_n,F,psip_n_fs_out)
        B_tor_out = F_out/R_fs_out
     
        # psi_pol_z selects ngrid_z points above and below Z0 and R_fs from psirz
        psi_pol_z = psirz[Z0_ind-ngrid_z:Z0_ind+ngrid_z,R_fs_ind_out+Rmag_ind]
        unif_z = np.linspace(Zgrid[Z0_ind-ngrid_z],Zgrid[Z0_ind+ngrid_z],10*ngrid_z)
        psi_pol_unifz = interp(Zgrid[Z0_ind-ngrid_z:Z0_ind+ngrid_z],psi_pol_z,unif_z)
        ###plt.plot(Zgrid[Z0_ind-ngrid_z:Z0_ind+ngrid_z],psi_pol_z,'x')
        ###plt.plot(unif_z,psi_pol_unifz,'.')
        ###plt.show()
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

        Z0_ind = Z0_ind+1

    R_fs_out_grid_r = np.flipud(R_fs_out_grid)
    B_pol_fs_out_r = np.flipud(B_pol_fs_out)
    B_tor_fs_out_r = np.flipud(B_tor_fs_out)

    B_tot_fs_out_r = np.sqrt(B_tor_fs_out_r**2+B_pol_fs_out_r**2)
    B_tot_fs_in = np.sqrt(B_tor_fs_in**2+B_pol_fs_in**2)

    theta_out_r = np.arctan(abs(Z_fs)/R_fs_out_grid_r)
    theta_in = -np.arctan(abs(Z_fs)/R_fs_in_grid)+np.pi

    R_fs = np.concatenate([R_fs_in_grid,R_fs_out_grid_r])
    B_pol = np.concatenate([B_pol_fs_in,B_pol_fs_out_r])
    B_tor = np.concatenate([B_tor_fs_in,B_tor_fs_out_r])
    B_tot = np.concatenate([B_tot_fs_in,B_tot_fs_out_r])
    theta_fs = np.concatenate([theta_in,theta_out_r])

    f=open('fs_'+efit_file_name+'.dat','w')
    f.write('# R B_pol B_tor B_tot theta\n')
    np.savetxt(f,np.column_stack((R_fs,B_pol,B_tor,B_tot,theta_fs)))
    f.close()

    return R_fs, B_pol, B_tor, B_tot
