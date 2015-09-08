from pylab import *
from sys import argv,exit,stdout
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
from finite_differences import *
from interp import *

def read_EFIT_file(efit_file_name):

    f = open(efit_file_name,'r')
    eqdsk=f.readlines()
    #print 'Header: %s' %eqdsk[0]
    #set resolutions
    nw=int(eqdsk[0].split()[-2])
    nh=int(eqdsk[0].split()[-1])
    #psi-width, number of flux surfaces around position of interest
    pw=(nw/8/2)*2 
    print 'EFIT file Resolution: %d x %d' %(nw,nh)

    entrylength=16
    #note: here rmin is rleft from EFIT
    try:
        rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
    except:
        entrylength=15
        try:
            rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
        except:
            exit('Error reading EQDSK file, please check format!')

    rmag,zmag,psiax,psisep,Bctr=[float(eqdsk[2][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[2])/entrylength)]
    dum,psiax2,dum,rmag2,dum=[float(eqdsk[3][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[3])/entrylength)]
    zmag2,dum,psisep2,dum,dum=[float(eqdsk[4][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[4])/entrylength)]
    if rmag!=rmag2: sys.exit('Inconsistent rmag: %7.4g, %7.4g' %(rmag,rmag2))
    if psiax2!=psiax: sys.exit('Inconsistent psiax: %7.4g, %7.4g' %(psiax,psiax2))
    if zmag!=zmag2: sys.exit('Inconsistent zmag: %7.4g, %7.4g' %(zmag,zmag2) )
    if psisep2!=psisep: sys.exit('Inconsistent psisep: %7.4g, %7.4g' %(psisep,psisep2))

    print "rmag", rmag
    print "zmag", zmag
    print "psiax", psiax
    print "psisep", psisep
    print "Bctr", Bctr
    # (R,Z) grid on which psi_pol is written
    Rgrid = np.arange(nw)/float(nw-1)*rdim+rmin
    print "rdim",rdim
    print "rmin",rmin
    print "first few Rgrid points", Rgrid[0:6]
    print "last few Rgrid points", Rgrid[-7:-1]
    Zgrid = np.arange(nh)/float(nh-1)*zdim+(zmid-zdim/2.0)
    print "zdim",zdim
    print "zmid",zmid
    print "first few Zgrid points", Zgrid[0:6]
    print "last few Zgrid points", Zgrid[-7:-1]

    # F, p, ffprime, pprime, q are written on uniform psi_pol grid
    # uniform grid of psi_pol~[psiax,psisep], resolution=nw
    F=empty(nw,dtype=float)
    p=empty(nw,dtype=float)
    ffprime=empty(nw,dtype=float)
    pprime=empty(nw,dtype=float)
    qpsi=empty(nw,dtype=float)
    # psi_pol is written on uniform (R,Z) grid (res=nw(R)*nh(Z))
    psirz_1d=empty(nw*nh,dtype=float)
    
    start_line=5
    lines=range(nw/5)
    if nw%5!=0: lines=range(nw/5+1)
    for i in lines:
        n_entries=len(eqdsk[i+start_line])/entrylength
        F[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    for i in lines:
        n_entries=len(eqdsk[i+start_line])/entrylength
        p[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    for i in lines:
        n_entries=len(eqdsk[i+start_line])/entrylength
        ffprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    for i in lines:
        n_entries=len(eqdsk[i+start_line])/entrylength
        pprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    lines_twod=range(nw*nh/5)
    if nw*nh%5!=0: lines_twod=range(nw*nh/5+1)
    for i in lines_twod:
        n_entries=len(eqdsk[i+start_line])/entrylength
        psirz_1d[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1
    psirz=psirz_1d.reshape(nh,nw)

    for i in lines:
        n_entries=len(eqdsk[i+start_line])/entrylength
        qpsi[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    # even grid of psi_pol, on which all 1D fields are defined
    psi_pol = linspace(psiax,psisep,nw)
    # return data read from efit file 
    # psi_pol grid, (R,Z) grid, F, p, ffprime, pprime, psirz, qpsi
    return psi_pol, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw

def calc_rho_tor(psi_pol, qpsi, nw):

    #create rho_tor grid on even psi_pol grid
    interpol_order = 3 
    q_spl_psi = US(psi_pol, qpsi, k=interpol_order, s=1e-5)
    psi_pol_fine = linspace(psi_pol[0], psi_pol[-1], nw*10)
    psi_tor_fine = empty((nw*10),dtype=float)
    psi_tor_fine[0] = 0.
    for i in range(1,nw*10):
        x=psi_pol_fine[:i+1]
        y=q_spl_psi(x)
        psi_tor_fine[i]=trapz(y,x)

    rho_tor_fine=sqrt(psi_tor_fine/psi_tor_fine[-1])
    rho_tor_spl=US(psi_pol_fine, rho_tor_fine, k=interpol_order, s=1e-5)
    # rho_tor grid (not uniform, on even grid of psi_pol) of resolution=nw 
    rho_tor=rho_tor_spl(psi_pol)

    return rho_tor_spl, psi_pol, rho_tor

def calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psi_pol, F):
    
    # Z0_ind is the index of Zgrid of midplane
    Z0_ind = np.argmin(np.abs(Zgrid-zmag))
    # psi_midplane is psi_pol at midplane on even Rgrid
    psi_pol_mp = psirz[Z0_ind,:]
    # Rmag_ind is index of Rgrid at rmag
    Rmag_ind = np.argmin(np.abs(Rgrid - rmag))
    print "rmag",rmag
    print "Rmag_ind",Rmag_ind
    print "Rgrid[Rmag_ind]",Rgrid[Rmag_ind]
    print "psi_pol_mp[Rmag_ind]",psi_pol_mp[Rmag_ind]
    print "psi_pol[-1]", psi_pol[-1]
    temp = np.copy(psi_pol_mp)
    temp[0:Rmag_ind] = 0.
    psi_pol_obmp = psi_pol_mp[Rmag_ind:]
    psisep_ind = np.argmin(np.abs(psi_pol_obmp-psi_pol[-1]))+7
    print "psisep_ind", psisep_ind
    print "psi_pol_obmp[psisep_ind]", psi_pol_obmp[psisep_ind]
    psi_pol_obmp = psi_pol_obmp[:psisep_ind]
    # B_pol is d psi_pol/ d R * (1/R)
    B_pol = fd_d1_o4(psi_pol_obmp, Rgrid[Rmag_ind:Rmag_ind+psisep_ind])/Rgrid[Rmag_ind:Rmag_ind+psisep_ind]
    # convert F(on even psi_pol grid) to F(on even R grid)
    F_evenR = interp(psi_pol, F, psi_pol_obmp) 
    # B_tor = F/R
    B_tor = F_evenR/Rgrid[Rmag_ind:Rmag_ind+psisep_ind]

    return psi_pol_obmp, B_pol, B_tor
