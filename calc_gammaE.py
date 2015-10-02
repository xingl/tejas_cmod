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
    v_d = Er/B_tot
    omega_i = e*B_tot/M
    omega_star = v_d*0.3*omega_i/v_thi
    omega_star_n = omega_star/(v_thi/a)
    gammaE = fd_d1_o4(Er/B_tot,R)/(v_thi/a)

    return omega_star_n,gammaE

def calc_gammaE_niti(R,ni,ti,ne,te,B_tor,B_pol,a):
    
    # assume n,p,B_tor,B_pol,Er,rho_tor are on uniform grid of R
    e = 1.6*10**(-19)
    M = 3.3*10**(-27)
    #a = 2.6863997038399379E-01
    p_i = ni*ti#+ne*te
    Er = fd_d1_o4(p_i,R)/ni/e
    B_tot = np.sqrt(B_tor**2+B_pol**2)
    v_thi = np.sqrt(ti/M)
    gammaE = fd_d1_o4(Er/B_tot,R)/(v_thi/a)
    p_tot = ni*ti+ne*te
    Er = fd_d1_o4(p_tot,R)/ni/e
    v_d = Er/B_tot
    omega_i = e*B_tot/M
    omega_star = v_d*0.3*omega_i/v_thi
    omega_star_n = omega_star/(v_thi/a)

    return omega_star_n,gammaE

def calc_gammaE_Er(R,ti,te,B_tor,B_pol,Er,a):

    # assume ti, B_tor, B_pol, Er are on uniform grid of R
    M = 3.3*10**(-27)
    #a = 2.6863997038399379E-01
    B_tot = np.sqrt(B_tor**2+B_pol**2)
    v_thi = np.sqrt(ti/M)
    v_E = Er/B_tot
    gammaE = fd_d1_o4(Er/B_tot,R)/(v_thi/a)
    cs = np.sqrt(te/M)
    Mach = Er/B_tot/cs
       
    return v_E,Mach,gammaE

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

def calc_gammaE_gene(R,rhot,te,B_pol,Er,q,a):

    # assume ti, B_tor, B_pol, Er are on uniform grid of R
    M = 3.3*10**(-27)
    uni_rhot = np.linspace(rhot[0],rhot[-1],1000)
    te_unirhot = interp(rhot,te,uni_rhot)
    Bp_unirhot = interp(rhot,B_pol,uni_rhot)
    Er_unirhot = interp(rhot,Er,uni_rhot)
    R_unirhot = interp(rhot,R,uni_rhot)
    q_unirhot = interp(rhot,q,uni_rhot)
    omega_tor = Er_unirhot/(R_unirhot*Bp_unirhot)
    cs = np.sqrt(te_unirhot/M)
    gammaE = (uni_rhot/q_unirhot)*fd_d1_o4(omega_tor,uni_rhot)/(cs/a)

    return uni_rhot, gammaE, omega_tor
