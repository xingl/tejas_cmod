from finite_differences import *
from interp import *
import matplotlib.pyplot as plt
import numpy as np

def calc_gammaE_p(R,ni,ti,p,B_tor,B_pol,a):

    # assume p, B_tor, B_pol, ni, ti are on uniform grid of R
    e = 1.6*10**(-19)
    M = 3.3*10**(-27)
    #a = 2.6863997038399379E-01
    Er = fd_d1_o4(p,R)/ni/e
    B_tot = np.sqrt(B_tor**2+B_pol**2)
    v_thi = np.sqrt(ti/M)
    gammaE = fd_d1_o4(Er/B_tot,R)/(v_thi/a)

    return gammaE

def calc_gammaE_niti(R,ni,ti,ne,te,B_tor,B_pol,a):
    
    # assume n,p,B_tor,B_pol,Er,rho_tor are on uniform grid of R
    e = 1.6*10**(-19)
    M = 3.3*10**(-27)
    #a = 2.6863997038399379E-01
    p = ni*ti#+ne*te
    Er = fd_d1_o4(p,R)/ni/e
    B_tot = np.sqrt(B_tor**2+B_pol**2)
    v_thi = np.sqrt(ti/M)
    gammaE = fd_d1_o4(Er/B_tot,R)/(v_thi/a)

    return gammaE

def calc_gammaE_Er(R,ti,B_tor,B_pol,Er,a):

    # assume ti, B_tor, B_pol, Er are on uniform grid of R
    M = 3.3*10**(-27)
    #a = 2.6863997038399379E-01
    B_tot = np.sqrt(B_tor**2+B_pol**2)
    v_thi = np.sqrt(ti/M)
    gammaE = fd_d1_o4(Er/B_tot,R)/(v_thi/a)
    
    return gammaE

def calc_gammaE_hb(R,psi_pol,ti,B_tor,B_pol,Er,rho_tor,a):
 
    # assume R,Er,B_tor,B_pol,ti,rho_tor are on psi_pol(not normalized)
    M = 3.3*10**(-27)
    #a = 2.6863997038399379E-01
    uni_psip = np.linspace(psi_pol[0],psi_pol[-1],500)
    R_unipsip = interp(psi_pol,R,uni_psip)
    Er_unipsip = interp(psi_pol,Er,uni_psip)
    Bt_unipsip = interp(psi_pol,B_tor,uni_psip)
    Bp_unipsip = interp(psi_pol,B_pol,uni_psip)
    ti_unipsip = interp(psi_pol,ti,uni_psip)
    v_thi = np.sqrt(ti_unipsip/M)
    B_tot = np.sqrt(Bt_unipsip**2+Bp_unipsip**2)
    omega_t = Er_unipsip/R_unipsip/Bp_unipsip
    gammaE = (R_unipsip*Bp_unipsip)**2*fd_d1_o4(omega_t,uni_psip)/B_tot
    gammaE = gammaE/(v_thi/a)
    rhot_unipsip = interp(psi_pol,rho_tor,uni_psip)

    return gammaE, omega_t, rhot_unipsip

