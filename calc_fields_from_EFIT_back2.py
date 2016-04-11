from scipy import interpolate
from read_EFIT_file import *
from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from finite_differences import *
from interp import *

#efit_file_name = 'efit_base'
#efit_file_name = 'efit_freebdry'
#efit_file_name = 'efit_ITER2C_I154_prof16_sym'
efit_file_name = 'efit_CRBS_ITER_n56_3a_freebdry'
flux_surface = '0.97'

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)
Z0_ind = np.argmin(abs(Zgrid-zmag))
print Z0_ind
psirz_lower = psirz[:Z0_ind+1,:].copy()
print 'shape of psirz_lower'
print np.shape(psirz_lower)
fs=float(flux_surface)
nR = len(Rgrid)
Zgrid_lower = Zgrid[0:Z0_ind+1]
nZ = len(Zgrid_lower)

print 'psiax = ', psiax
print 'psisep = ', psisep
print 'len(Rgrid)', 'len(Zgrid)'
print nR, nZ
print 'rmag','zmag'
print rmag,zmag
print Zgrid[Z0_ind]
print 'shape of psirz'
print np.shape(psirz)

psirz_spl = interpolate.RectBivariateSpline(Zgrid_lower,Rgrid,psirz_lower)
#newpsirz = psirz_spl(Zgrid,Rgrid)
#diff_psirz = psirz-newpsirz
#print max(abs(diff_psirz.flatten()))

Bp_R_grid = np.empty(np.shape(psirz_lower))
for i in range(nZ):
    #print i
    Bp_R_grid_this=fd_d1_o4(psirz_lower[i,:].flatten(),Rgrid)/Rgrid
    Bp_R_grid[i,:]=Bp_R_grid_this.copy()
    Bp_R_grid[i,0]=(psirz_lower[i,0]-psirz_lower[i,1])/(Rgrid[0]-Rgrid[1])
    Bp_R_grid[i,1]=(psirz_lower[i,1]-psirz_lower[i,2])/(Rgrid[1]-Rgrid[2])
    Bp_R_grid[i,nR-1]=(psirz_lower[i,nR-1]-psirz_lower[i,nR-2])/(Rgrid[nR-1]-Rgrid[nR-2])
    Bp_R_grid[i,nR-2]=(psirz_lower[i,nR-2]-psirz_lower[i,nR-3])/(Rgrid[nR-2]-Rgrid[nR-3])
    #print Bp_R_grid[i,nR-3:nR]
    #print Bp_R_grid[i,0:3]
Bp_Z_grid = np.empty(np.shape(psirz_lower))
for i in range(nR):
    #print i
    Bp_Z_grid_this=fd_d1_o4(psirz_lower[:,i].flatten(),Zgrid_lower)/Rgrid[i]
    Bp_Z_grid[:,i]=Bp_Z_grid_this.copy()
    Bp_Z_grid[0,i]=(psirz_lower[0,i]-psirz_lower[1,i])/(Zgrid_lower[0]-Zgrid_lower[1])
    Bp_Z_grid[1,i]=(psirz_lower[1,i]-psirz_lower[2,i])/(Zgrid_lower[1]-Zgrid_lower[2])
    Bp_Z_grid[nZ-1,i]=(psirz_lower[nZ-1,i]-psirz_lower[nZ-2,i])/(Zgrid_lower[nZ-1]-Zgrid_lower[nZ-2])
    Bp_Z_grid[nZ-2,i]=(psirz_lower[nZ-2,i]-psirz_lower[nZ-3,i])/(Zgrid_lower[nZ-2]-Zgrid_lower[nZ-3])
    #print Bp_Z_grid[nZ-3:nZ,i]
    #print Bp_Z_grid[0:3,i]

Bp_R_spl = interpolate.RectBivariateSpline(Zgrid_lower,Rgrid,Bp_R_grid)
Bp_Z_spl = interpolate.RectBivariateSpline(Zgrid_lower,Rgrid,Bp_Z_grid)

Bp_tot_grid = np.sqrt(Bp_R_grid**2+Bp_Z_grid**2)
Bp_obmp = Bp_tot_grid[Z0_ind,:]
#Bp_obmp = np.sqrt(Bp_R_grid[Z0_ind,:].flatten()**2+Bp_Z_grid[Z0_ind,:].flatten()**2)
f = open('Bp_obmp.dat','w')
f.write('# R Bp_t\n # \n')
np.savetxt(f,np.column_stack((Rgrid,Bp_obmp)))
f.close()
#rbs = np.genfromtxt('rbsProfs')
#plt.plot(rbs[:,24],rbs[:,25],'r',label='rbsProfs')
#plt.plot(Rgrid,Bp_obmp,'.',label='spline')
##plt.plot(Rgrid,Bp_R_grid[Z0_ind,:],label='spline')
#plt.legend(loc=2)
#plt.title('outboard midplane')
#plt.ylabel('B_pol')
#plt.xlabel('R')
#plt.title(efit_file_name)
#plt.show()

psip_mp = psirz_lower[Z0_ind,:].copy()
R0_ind = np.argmin(psip_mp)
print 'R0_ind = ', R0_ind
psip_obmp = psip_mp[R0_ind:].copy()
R_fs = np.argmin(abs(psip_obmp-(fs*(psisep-psiax)+psiax)))
print 'R_fs = ', R_fs
Bp_tan = Bp_tot_grid[Z0_ind-50:Z0_ind,R_fs+R0_ind].copy()
plt.plot(Zgrid_lower[Z0_ind-50:Z0_ind],Bp_tan.flatten(),'.')
plt.ylabel('B_pol')
plt.xlabel('Z')
plt.title(efit_file_name)
plt.show()

plt.plot(Rgrid,Bp_tot_grid[Z0_ind-2,:],label='Zmp-2')
plt.plot(Rgrid,Bp_tot_grid[Z0_ind-1,:],label='Zmp-1')
plt.plot(Rgrid,Bp_tot_grid[Z0_ind,:],label='Zmp')
plt.ylabel('B_pol')
plt.xlabel('R')
plt.legend(loc=2)
plt.title(efit_file_name)
plt.show()

tol = 1.0E-06
l_1 = 1.8
l_2 = 1.7
ntheta = 1000
l_grid = np.empty(ntheta)
R_fs = np.empty(ntheta)
Z_fs = np.empty(ntheta)
Bp_R = np.empty(ntheta)
Bp_Z = np.empty(ntheta)
Bp = np.empty(ntheta)

for i in range(ntheta):
    theta = np.pi*i/ntheta
    #print i, theta
    R_1 = rmag+l_1*np.cos(theta)
    Z_1 = zmag-l_1*np.sin(theta)
    psi_1 = (psirz_spl(Z_1,R_1)[0][0]-psiax)/(psisep-psiax) - fs
    #print R_1,Z_1,psi_1
    R_2 = rmag+l_2*np.cos(theta)
    Z_2 = zmag-l_2*np.sin(theta)
    psi_2 = (psirz_spl(Z_2,R_2)[0][0]-psiax)/(psisep-psiax) - fs
    #print R_2,Z_2,psi_2
    istep = 0
    while abs(psi_1)>tol and abs(psi_2)>tol and istep<20:
        this_l = l_1-psi_1*(l_1-l_2)/(psi_1-psi_2)
	this_R = rmag+this_l*np.cos(theta)
	this_Z = zmag-this_l*np.sin(theta)
	this_psi = (psirz_spl(this_Z,this_R)[0][0]-psiax)/(psisep-psiax) - fs
	psi_2 = psi_1
	psi_1 = this_psi
	l_2 = l_1
	l_1 = this_l
	istep = istep +1
	#print istep
	#print this_l, this_psi 
    l_grid[i] = l_1
    R_fs[i] = rmag+l_grid[i]*np.cos(theta)
    Z_fs[i] = zmag-l_grid[i]*np.sin(theta)
    Bp_R[i] = Bp_R_spl(Z_fs[i],R_fs[i])[0][0]
    Bp_Z[i] = Bp_Z_spl(Z_fs[i],R_fs[i])[0][0]
    Bp[i] = np.sqrt(Bp_R[i]**2+Bp_Z[i]**2)
#plt.plot(R_fs,Z_fs)
#plt.show()

#data = np.genfromtxt('Bp.dat')
#plt.plot(data[:,0],data[:,1])
#plt.plot(R_fs,Bp)
#plt.show()

F_spl = interpolate.UnivariateSpline(psip_n,F)
F_fs = F_spl(fs)
print fs,F_fs
Bt = F_fs/R_fs
#plt.plot(data[:,0],data[:,2])
#plt.plot(R_fs,Bt)
#plt.show()

f = open('fs.dat','w')
f.write('# R Z Bp_R Bp_Z\n # \n')
np.savetxt(f,np.column_stack((R_fs,Z_fs,Bp_R,Bp_Z)))
f.close()
