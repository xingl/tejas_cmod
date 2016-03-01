from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *
from read_EFIT_file import *
from fields_along_fluxsurface import *
from read_rbsProfs import *

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = args[0]
rbsProfs_file_name = args[1]
fs = args[2]

R_fs, B_pol, B_tor, B_tot = fields_along_fluxsurface(efit_file_name,fs)

plt.plot(R_fs,B_pol,label='B_pol')
plt.plot(R_fs,abs(B_tor),label='B_tor')
plt.plot(R_fs,B_tot,label='B_tot')
plt.title(efit_file_name+' at flux surface of '+fs+' +- 0.01')
plt.xlabel('R')
plt.ylabel('B')
plt.legend()
plt.show()

R_grid,tprime_fs,fprime_fs,gamE_n_fs = read_rbsProfs(rbsProfs_file_name,fs,efit_file_name)

plt.plot(R_grid,tprime_fs,label='tprime')
plt.plot(R_grid,fprime_fs,label='fprime')
plt.xlabel('R')
plt.legend(loc=2)
plt.title(rbsProfs_file_name+' fs= '+str(fs))
plt.show()
plt.plot(R_grid,gamE_n_fs,label='gammaE')
plt.xlabel('R')
plt.legend(loc=4)
plt.title(rbsProfs_file_name+' fs= '+str(fs))
plt.show()

