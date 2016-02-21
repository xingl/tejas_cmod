import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences import *

rbsProfs_data = np.genfromtxt('rbsProfs')
fs = 0.98
a = rbsProfs_data[-1,22]

rhot_data = rbsProfs_data[:,0]
fs_ind = np.argmin(abs(rhot_data-fs))
rhot = rhot_data[fs_ind-10:fs_ind+10]
psip = rbsProfs_data[fs_ind-10:fs_ind+10,1]
n = rbsProfs_data[fs_ind-10:fs_ind+10,3]
Ti = rbsProfs_data[fs_ind-10:fs_ind+10,4]
Te = rbsProfs_data[fs_ind-10:fs_ind+10,5]
R = rbsProfs_data[fs_ind-10:fs_ind+10,24]
B_pol = rbsProfs_data[fs_ind-10:fs_ind+10,25]
B_tor = rbsProfs_data[fs_ind-10:fs_ind+10,26]
gamE = rbsProfs_data[fs_ind-10:fs_ind+10,9]
#print R[10],B_pol[10],B_tor[10]
unif_R = np.linspace(R[0],R[-1],100)
Ti_unifR = interp(R,Ti,unif_R)
tprime = a*fd_d1_o4(Ti_unifR,unif_R)/Ti_unifR
n_unifR = interp(R,n,unif_R)
fprime = a*fd_d1_o4(n_unifR,unif_R)/n_unifR
gamE_unifR = interp(R,gamE,unif_R)

R_fs_ind = np.argmin(abs(unif_R-R[10]))
#print tprime[R_fs_ind],fprime[R_fs_ind],gamE_prime[R_fs_ind]

tprime_obmp = abs(tprime[R_fs_ind])
fprime_obmp = abs(fprime[R_fs_ind])
gamE_obmp = abs(gamE_unifR[R_fs_ind])
R_obmp = R[10]
B_pol_obmp = abs(B_pol[10])
B_tor_obmp = abs(B_tor[10])
B_tot_obmp = np.sqrt(B_pol_obmp**2+B_tor_obmp**2)

fs_data = np.genfromtxt('fs.dat')
R_grid = fs_data[:,0]
B_pol = abs(fs_data[:,1])
B_tor = abs(fs_data[:,2])
B_tot = fs_data[:,3]

tprime_fs = R_grid*B_pol*tprime_obmp/R_obmp/B_pol_obmp
fprime_fs = R_grid*B_pol*fprime_obmp/R_obmp/B_pol_obmp
gammaE_fs = (R_grid*B_pol)**2/B_tot*gamE_obmp*B_tot_obmp/(R_obmp*B_pol_obmp)**2
plt.plot(R_grid,tprime_fs,'.',label='tprime')
plt.plot(R_grid,fprime_fs,'.',label='fprime')
plt.xlabel('R')
plt.legend()
plt.show()
plt.plot(R_grid,gammaE_fs,'.')
plt.ylabel('gammaE')
plt.xlabel('R')
plt.show()
