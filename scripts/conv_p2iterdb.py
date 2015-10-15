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
from read_iterdb_file import *
from w_iterdb import *

e = 1.6*10**(-19)
a = 2.6863997038399379E-01

parser = op.OptionParser()
parser.add_option('--shift_Ti','-s',action='store_const',const=1,help='shift Ti profile out by 0.005',default=0)
parser.add_option('--output_iterdb','-o',action='store_const',const=1,help='write out iterdb file from p file',default=0)
options,args = parser.parse_args()
if len(args) == 4:
    efit_file_name = args[0]
    p_file_name = args[1]
    Er_file_name = args[2]
    iterdb_filename = args[3]
elif len(args) == 3:
    efit_file_name = args[0]
    p_file_name = args[1]
    Er_file_name = args[2]
else:
    exit("""
Please include EFIT, p, Er, iterdb(optional) files as arguments (in that order).
     """)
shift_Ti = options.shift_Ti
output_iterdb =options.output_iterdb

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file(efit_file_name)

rho_tor_spl, rhot_n = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw)

psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F,nw,psip_n)

psip_obmp = psip_n_obmp*(psisep-psiax)+psiax
rhot_n_obmp = rho_tor_spl(psip_obmp)

### read from p-file
### ne, ni are in the unit of 10^20 m^(-3)
### te, ti are in the unit of KeV
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

rhotp_obmp = interp(psi0,rhot0,psip_n_obmp)
q_obmp = interp(psip_n, qpsi, psip_n_obmp)

if shift_Ti:
    shift_Er = True
else:
    shift_Er = False
Er_shift = 0.005  #Shift in poloidal flux coord
psi0_Er, Er, Er_error = read_Er(Er_file_name,shift_Er,Er_shift)
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

rhot_gene,gammaE_gene,omega_tor_gene = calc_gammaE_gene(R_obmp[psi_Er_f:psi_Er_l+1], rhot_n_obmp[psi_Er_f:psi_Er_l+1], te_obmp[psi_Er_f:psi_Er_l+1], B_pol[psi_Er_f:psi_Er_l+1], Er_obmp, q_obmp[psi_Er_f:psi_Er_l+1], a)

omega_tor_full = interp(rhot_gene,omega_tor_gene,rhot0)

if len(args) == 4:
    ### read from iterdb file
    ### te, ti are in the unit of ev
    ### ne, ni are in the unit of m^(-3)
    rhot_idb, te_idb, ti_idb, ne_idb, ni_idb, nb_idb, vrot_idb = read_iterdb_file(iterdb_filename)
    psip_idb = interp(rhot_n,psip_n,rhot_idb)
    psipne, ne, psipte, te, psipni, ni, psipti, ti = read_cmod_pfile_raw(p_file_name)
    
    plt.plot(psip_idb,te_idb,'x',label='te iterdb')
    plt.plot(psipte,te*1E03,'r.',label='te p-file')
    plt.xlabel('psip_n')
    plt.legend()
    plt.show()

    plt.plot(psip_idb,ti_idb,'x',label='ti iterdb')
    plt.plot(psipti,ti*1E03,'r.',label='ti p-file')
    plt.xlabel('psip_n')
    plt.legend()
    plt.show()

    plt.plot(psip_idb,ne_idb,'x',label='ne iterdb')
    plt.plot(psipne,ne*1E20,'r.',label='ne p-file')
    plt.xlabel('psip_n')
    plt.legend()
    plt.show()

    plt.plot(psip_idb,ni_idb,'x',label='ni iterdb')
    plt.plot(psipni,ni*1E20,'r.',label='ni p-file')
    plt.xlabel('psip_n')
    plt.legend()
    plt.show()

    plt.plot(rhot0,omega_tor_full,'x',label='p-file')
    plt.plot(rhot_idb,vrot_idb,'.',label='iterdb')
    plt.axis((0.9,1.0,-150000,20000))
    plt.legend()
    plt.show()

if output_iterdb:
    file_out_base = 'profiles_wb' 
    base_number = '1120907032.01012'
    output_iterdb(rhot0,ne,te/e,ni,ti/e,file_out_base+base_number,base_number,'9999',vrot=omega_tor_full,nimp=nz)

    f = open(file_out_base+'_gene_e.dat','w')
    f.write('# rhot_n psip_n Te(eV) ne(m^-3) \n # \n')
    np.savetxt(f,np.column_stack((rhot0,psi0,te/e,ne)))
    f.close()

    f = open(file_out_base+'_gene_i.dat','w')
    f.write('# rhot_n psip_n Ti(eV) ni(m^-3) \n # \n')
    np.savetxt(f,np.column_stack((rhot0,psi0,ti/e,ni)))
    f.close()
