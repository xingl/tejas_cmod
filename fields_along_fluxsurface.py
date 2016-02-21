from pylab import *
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *
from read_EFIT_file import *

efit_file_name = 'g_new_901_901_1415'
fs = 0.98

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)


Z0_ind = np.argmin(np.abs(Zgrid-zmag))
Rmag_ind = np.argmin(np.abs(Rgrid - rmag))

psi_pol_mp = psirz[Z0_ind,:]
psi_pol_out = psi_pol_mp[Rmag_ind:].copy()
unif_R_out_1 = np.linspace(Rgrid[Rmag_ind],Rgrid[-1],nw*10)
psip_n_temp = (psi_pol_out-psiax)/(psisep-psiax)
psip_n_mp_out = interp(Rgrid[Rmag_ind:],psip_n_temp,unif_R_out_1)
F_mp_out = interp(psip_n,F,psip_n_mp_out)
B_tor_mp_out=F_mp_out/unif_R_out_1

psi_pol_in = psi_pol_mp[:Rmag_ind].copy()
unif_R_in_1 = np.linspace(Rgrid[0],Rgrid[Rmag_ind],nw*10)
psip_n_temp = (psi_pol_in-psiax)/(psisep-psiax)
psip_n_mp_in = interp(Rgrid[:Rmag_ind],psip_n_temp,unif_R_in_1)
F_mp_in = interp(psip_n,F,psip_n_mp_in)
B_tor_mp_in=F_mp_in/unif_R_in_1
#psip_n_1,R_1,B_p,B_t=calc_B_fields(Rgrid,rmag,Zgrid,zmag,psirz,psiax,psisep,F,nw,psip_n)

interpol_order = 3
Bt_mp_out = US(unif_R_out_1,B_tor_mp_out,k=interpol_order,s=1e-5)
Bt_mp_in = US(unif_R_in_1,B_tor_mp_in,k=interpol_order,s=1e-5)
#plt.plot(unif_R_out_1,B_tor_mp_out)
#plt.plot(unif_R_in_1,B_tor_mp_in)
#plt.show()

Z_fs = np.empty(0) 
R_fs_out_grid = np.empty(0) 
B_pol_fs_out = np.empty(0)
R_fs_in_grid = np.empty(0)
B_pol_fs_in = np.empty(0)
theta_out = np.empty(0)
theta_in = np.empty(0)

while (np.min((psirz[Z0_ind,:]-psiax)/(psisep-psiax))<fs):

    psi_pol = psirz[Z0_ind,:]
    Rmag_ind = np.argmin(abs(psi_pol))
    psi_pol_out = psi_pol[Rmag_ind:].copy()
    psip_n_temp = (psi_pol_out-psiax)/(psisep-psiax)
    R_fs_ind_out = np.argmin(abs(psip_n_temp-fs))
    unif_R_out = np.linspace(Rgrid[Rmag_ind],Rgrid[-1],nw*10)
    psip_n_unifR_out = interp(Rgrid[Rmag_ind:],psip_n_temp,unif_R_out)
    psifs_ind_out = np.argmin(abs(psip_n_unifR_out-fs))
    psip_n_fs_out = psip_n_unifR_out[psifs_ind_out-3:psifs_ind_out+3].copy()
    R_fs_out = unif_R_out[psifs_ind_out-3:psifs_ind_out+3].copy()
    B_pol_Z_out = fd_d1_o4(psip_n_fs_out*(psisep-psiax)+psiax,R_fs_out)/R_fs_out

    psi_pol_z = psirz[Z0_ind-10:Z0_ind+10,R_fs_ind_out+Rmag_ind]
    unif_z = np.linspace(Zgrid[Z0_ind-10],Zgrid[Z0_ind+10],100)
    psi_pol_unifz = interp(Zgrid[Z0_ind-10:Z0_ind+10],psi_pol_z,unif_z)
    B_pol_R_out = fd_d1_o4(psi_pol_unifz,unif_z)/R_fs_out[3]
    z_fs_ind=np.argmin(abs(unif_z-Zgrid[Z0_ind]))
    B_pol_out = np.sqrt(B_pol_Z_out[3]**2+B_pol_R_out[z_fs_ind]**2) 
    
    psi_pol_in = psi_pol[:Rmag_ind].copy()
    psip_n_temp = (psi_pol_in-psiax)/(psisep-psiax)
    R_fs_ind_in = np.argmin(abs(psip_n_temp-fs))
    unif_R_in = np.linspace(Rgrid[0],Rgrid[Rmag_ind],nw*10)
    psip_n_unifR_in = interp(Rgrid[:Rmag_ind],psip_n_temp,unif_R_in)
    psifs_ind_in = np.argmin(abs(psip_n_unifR_in-fs))
    psip_n_fs_in = psip_n_unifR_in[psifs_ind_in-3:psifs_ind_in+3].copy()
    R_fs_in = unif_R_in[psifs_ind_in-3:psifs_ind_in+3].copy()
    B_pol_Z_in = fd_d1_o4(psip_n_fs_in*(psisep-psiax)+psiax,R_fs_in)/R_fs_in

    psi_pol_z = psirz[Z0_ind-10:Z0_ind+10,R_fs_ind_in]
    psi_pol_unifz = interp(Zgrid[Z0_ind-10:Z0_ind+10],psi_pol_z,unif_z)
    B_pol_R_in = fd_d1_o4(psi_pol_unifz,unif_z)/R_fs_in[3]
    B_pol_in = np.sqrt(B_pol_Z_in[3]**2+B_pol_R_in[z_fs_ind]**2) 
#    print B_pol_Z_out[3],B_pol_R_out[z_fs_ind],B_pol_out,R_fs_out[3]
#    print B_pol_Z_in[3],B_pol_R_in[z_fs_ind],B_pol_in,R_fs_in[3]
#    print Z0_ind

    Z_fs = np.append(Z_fs,Zgrid[Z0_ind])
    R_fs_out_grid = np.append(R_fs_out_grid,R_fs_out[3])
    B_pol_fs_out = np.append(B_pol_fs_out,B_pol_out)
    #B_pol_fs_out = np.append(B_pol_fs_out,B_pol_Z_out[3])
    #B_pol_fs_out = np.append(B_pol_fs_out,B_pol_R_out[z_fs_ind])
    R_fs_in_grid = np.append(R_fs_in_grid,R_fs_in[3])
    B_pol_fs_in = np.append(B_pol_fs_in,B_pol_in)
    #B_pol_fs_in = np.append(B_pol_fs_in,B_pol_Z_in[3])
    #B_pol_fs_in = np.append(B_pol_fs_in,B_pol_R_in[z_fs_ind])

    Z0_ind = Z0_ind-1

#    print B_pol_out[3],R_fs_out[3]
#    print B_pol_in[3],R_fs_in[3]
#    print Z0_ind


#plt.plot(R_fs_out_grid,B_pol_fs_out)
#plt.plot(R_fs_in_grid,B_pol_fs_in)
#plt.show()

B_tor_fs_out = Bt_mp_out(R_fs_out_grid)
B_tor_fs_in = Bt_mp_in(R_fs_in_grid)

#print R_fs_out_grid[0],B_pol_fs_out[0],B_tor_fs_out[0]

#plt.plot(Z_fs,abs(B_pol_fs_out),label='B_pol_out')
#plt.plot(Z_fs,abs(B_pol_fs_in),label='B_pol_in')
#plt.plot(Z_fs,abs(B_tor_fs_out),'.',label='B_tor_out')
#plt.plot(Z_fs,abs(B_tor_fs_in),'.',label='B_tor_in')
#plt.legend()
#plt.xlabel('Z')
#plt.ylabel('T')
#plt.show()

#B_tor_fs_out=interp(unif_R_out,B_tor_mp_out,R_fs_out_grid)
#B_tor_fs_in=interp(unif_R_in,B_tor_mp_in,R_fs_in_grid)
plt.plot(R_fs_out_grid,abs(B_pol_fs_out),label='B_pol_out')
plt.plot(R_fs_in_grid,abs(B_pol_fs_in),label='B_pol_in')
plt.plot(R_fs_out_grid,abs(B_tor_fs_out),'.',label='B_tor_out')
plt.plot(R_fs_in_grid,abs(B_tor_fs_in),'.',label='B_tor_in')
#plt.plot(unif_R_out_1,abs(B_tor_mp_out),label='B_tor_out_cmp')
#plt.plot(unif_R_in_1,abs(B_tor_mp_in),label='B_tor_in_cmp')
plt.legend()
plt.xlabel('R(m)')
plt.ylabel('B(T)')
plt.show()
B_tot_fs_out = np.sqrt(B_tor_fs_out**2+B_pol_fs_out**2)
B_tot_fs_in = np.sqrt(B_tor_fs_in**2+B_pol_fs_in**2)

R_fs = np.concatenate([R_fs_out_grid,R_fs_in_grid])
B_pol = np.concatenate([B_pol_fs_out,B_pol_fs_in])
B_tor = np.concatenate([B_tor_fs_out,B_tor_fs_in])
B_tot = np.concatenate([B_tot_fs_out,B_tot_fs_in])

f=open('fs.dat','w')
f.write('# R B_pol B_tor B_tot\n')
np.savetxt(f,np.column_stack((R_fs,B_pol,B_tor,B_tot)))
f.close
