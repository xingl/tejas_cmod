from scipy import interpolate
from read_EFIT_file import *
from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from finite_differences import *
from interp import *

#efit_file_name = 'g_new_901_901_1415'
#efit_file_name = 'FS2RZ_gEFIT'
#flux_surface = '0.97'
##
parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = args[0]
flux_surface = args[1]

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)
Z0_ind = np.argmin(abs(Zgrid-zmag))
print Z0_ind
#psirz_lower = psirz[:Z0_ind,:].copy()
psirz_lower = psirz.copy()
print 'shape of psirz_lower'
print np.shape(psirz_lower)
fs=float(flux_surface)
nR = len(Rgrid)
#Zgrid_lower = Zgrid[0:Z0_ind]
Zgrid_lower = Zgrid.copy()
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

Bp_R = np.empty(np.shape(psirz_lower))
for i in range(nZ):
    #print i
    Bp_R_this=fd_d1_o4(psirz_lower[i,:].flatten(),Rgrid)/Rgrid
    Bp_R[i,:]=Bp_R_this.copy()
    Bp_R[i,0]=(psirz_lower[i,0]-psirz_lower[i,1])/(Rgrid[0]-Rgrid[1])
    Bp_R[i,1]=(psirz_lower[i,1]-psirz_lower[i,2])/(Rgrid[1]-Rgrid[2])
    Bp_R[i,nR-1]=(psirz_lower[i,nR-1]-psirz_lower[i,nR-2])/(Rgrid[nR-1]-Rgrid[nR-2])
    Bp_R[i,nR-2]=(psirz_lower[i,nR-2]-psirz_lower[i,nR-3])/(Rgrid[nR-2]-Rgrid[nR-3])
    #print Bp_R[i,nR-3:nR]
    #print Bp_R[i,0:3]
Bp_Z = np.empty(np.shape(psirz_lower))
for i in range(nR):
    #print i
    Bp_Z_this=fd_d1_o4(psirz_lower[:,i].flatten(),Zgrid_lower)/Rgrid[i]
    Bp_Z[:,i]=Bp_Z_this.copy()
    Bp_Z[0,i]=(psirz_lower[0,i]-psirz_lower[1,i])/(Zgrid_lower[0]-Zgrid_lower[1])
    Bp_Z[1,i]=(psirz_lower[1,i]-psirz_lower[2,i])/(Zgrid_lower[1]-Zgrid_lower[2])
    Bp_Z[nZ-1,i]=(psirz_lower[nZ-1,i]-psirz_lower[nZ-2,i])/(Zgrid_lower[nZ-1]-Zgrid_lower[nZ-2])
    Bp_Z[nZ-2,i]=(psirz_lower[nZ-2,i]-psirz_lower[nZ-3,i])/(Zgrid_lower[nZ-2]-Zgrid_lower[nZ-3])
    #print Bp_Z[nZ-3:nZ,i]
    #print Bp_Z[0:3,i]
Bp_R_spl = interpolate.RectBivariateSpline(Zgrid_lower,Rgrid,Bp_R)
Bp_Z_spl = interpolate.RectBivariateSpline(Zgrid_lower,Rgrid,Bp_Z)

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
    theta = 2.*np.pi*i/ntheta
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
plt.plot(R_fs,Z_fs)
plt.title(efit_file_name)
plt.show()

#data = np.genfromtxt('Bp.dat')
data_base = np.genfromtxt('Bp_base_new.dat')
plt.plot(data_base[:,0],data_base[:,1],label='base')
#plt.plot(data[:,0],data[:,1],label='old')
plt.plot(R_fs,Bp,label='new')
plt.xlabel('R')
plt.ylabel('B_pol')
plt.title(efit_file_name)
plt.legend()
plt.show()

F_spl = interpolate.UnivariateSpline(psip_n,F)
F_fs = F_spl(fs)
print fs,F_fs
Bt = F_fs/R_fs
plt.plot(data_base[:,0],abs(data_base[:,2]),label='base')
#plt.plot(data[:,0],data[:,2],label='old')
plt.plot(R_fs,Bt,label='new')
plt.xlabel('R')
plt.ylabel('B_tor')
plt.title(efit_file_name)
plt.legend()
plt.show()

f = open('Bp_new.dat','w')
f.write('# R Bp Bt\n # \n')
np.savetxt(f,np.column_stack((R_fs,Bp,Bt)))
f.close()
