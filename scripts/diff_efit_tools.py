#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots: q(psi)
       P(psi)
       psi(R,Z)
       psi(R,Z) for psi = 0.9 and 1.0 along with grid points
-c: outputs rho_tor and rho_pol
-p: outputs R, psi, B_pol
-n: suppresses plots

Modified from extract_miller_from_eqdsk.py
"""

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate

parser=op.OptionParser(description='Tools for extracting information from an EFIT file.')
#parser.add_option('--rovera','-r',action='store_const',const=1)
parser.add_option('--conv','-c',action='store_const',const=1,help = 'Output rho_tor vs rho_pol')
parser.add_option('--binfo','-p',action='store_const',const=1,help = 'Ouptut R, psi_pol, B_pol, B_tor')
parser.add_option('--noplot','-n',action='store_const',const=1,help = 'Suppress plots')
parser.add_option('--contour','-r',action='store_const',const=1,help = 'Output contours') 
options,args=parser.parse_args()
write_rhotor_rhopol_conversion=options.conv
write_binfo=options.binfo
show_contour=options.contour
show_plots = not options.noplot
#use_r_a=options.rovera
if len(args) != 2:
    exit("""
Please include two efit file names as argument, lower resolution first."
    \n""")
    
filename_lowres=args[0]
file_lowres=open(filename_lowres,'r')

filename_highres=args[1]
file_highres=open(filename_highres,'r')


def find(val,arr):
    return argmin(abs(arr-val))

def fd_d1_o4(var,grid,mat=False):
    """Centered finite difference, first derivative, 4th order.
    var: quantity to be differentiated.
    grid: grid for var 
    mat: matrix for the finite-differencing operator. if mat=False then it is created"""

    if not mat:
        mat=get_mat_fd_d1_o4(len(var),grid[1]-grid[0])

    dvar=np.dot(mat,var)
    dvar[0]=0.0
    dvar[1]=0.0
    #dvar[2]=0.0
    dvar[-1]=0.0
    dvar[-2]=0.0
    #dvar[-3]=0.0
    return dvar 

def get_mat_fd_d1_o4(size,dx,plot_matrix=False):
    """Creates matrix for centered finite difference, first derivative, 4th order.
    size: size of (number of elements in) quantity to be differentiated
    dx: grid spacing (for constant grid)."""

    prefactor=1.0/(12.0*dx)
    mat=np.zeros((size,size),dtype='float')    
    for i in range(size):
        if i-1 >= 0:
            mat[i,i-1]=-8
        if i-2 >= 0:
            mat[i,i-2]=1
        if i+1 <= size-1:
            mat[i,i+1]=8
        if i+2 <= size-1:
            mat[i,i+2]=-1
   
    mat=prefactor*mat

    if plot_matrix:
        plt.contourf(mat,50)
        plt.colorbar()
        plt.show()

    return mat

def interp(xin,yin,xnew):
    """
    xin: x variable input
    yin: y variable input
    xnew: new x grid on which to interpolate
    yout: new y interpolated on xnew
    """

    #splrep returns a knots and coefficients for cubic spline
    rho_tck = interpolate.splrep(xin,yin)
    #Use these knots and coefficients to get new y
    yout = interpolate.splev(xnew,rho_tck,der=0)

    return yout

def full_interp(func_xin,xin,xconv,yconv,yout):
    """
    Takes function func_xin on grid xin and outputs the function on yout grid
    func_xin: function to interpolate
    xin: grid corresponding to func_xin
    xconv: xgrid for conversion
    yconv: ygrid for conversion
    yout: output grid
    """

    #If necessary truncate func_xin onto correct range
    if xin[0] < xconv[0]:
        low_index = np.argmin(abs(xconv-xin[0]))
    else:
        low_index = 0
    if xin[-1] > xconv[-1]:
        high_index = np.argmin(abs(xconv-xin[-1]))
    else:
        high_index = -1

    if high_index == -1:
        func_xin = func_xin[low_index:]
    else:
        func_xin = func_xin[low_index:high_index]

    func_xconv = interp(xin,func_xin,xconv)
    func_yout = interp(yconv,func_xconv,yout)
    
    return func_yout
    
eqdsk_lowres=file_lowres.readlines()
print 'Header: %s' %eqdsk_lowres[0]
#set resolutions
nw_lowres=int(eqdsk_lowres[0].split()[-2]);nh_lowres=int(eqdsk_lowres[0].split()[-1])
pw_lowres=(nw_lowres/8/2)*2 #psi-width, number of flux surfaces around position of interest
print 'Resolution: %d x %d' %(nw_lowres,nh_lowres)

entrylength=16
#note: here rmin_lowres is rleft from EFIT
try:
    rdim_lowres,zdim_lowres,rctr_lowres,rmin_lowres,zmid_lowres=[float(eqdsk_lowres[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_lowres[1])/entrylength)]
except:
    entrylength=15
    try:
        rdim_lowres,zdim_lowres,rctr_lowres,rmin_lowres,zmid_lowres=[float(eqdsk_lowres[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_lowres[1])/entrylength)]
    except:
        exit('Error reading EQDSK file, please check format!')

rmag_lowres,zmag_lowres,psiax_lowres,psisep_lowres,Bctr_lowres=[float(eqdsk_lowres[2][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_lowres[2])/entrylength)]
dum_lowres,psiax_lowres2,dum_lowres,rmag_lowres2,dum_lowres=[float(eqdsk_lowres[3][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_lowres[3])/entrylength)]
zmag_lowres2,dum_lowres,psisep_lowres2,dum_lowres,dum_lowres=[float(eqdsk_lowres[4][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_lowres[4])/entrylength)]
if rmag_lowres!=rmag_lowres2: sys.exit('Inconsistent rmag_lowres: %7.4g, %7.4g' %(rmag_lowres,rmag_lowres2))
if psiax_lowres2!=psiax_lowres: sys.exit('Inconsistent psiax_lowres: %7.4g, %7.4g' %(psiax_lowres,psiax_lowres2))
if zmag_lowres!=zmag_lowres2: sys.exit('Inconsistent zmag_lowres: %7.4g, %7.4g' %(zmag_lowres,zmag_lowres2) )
if psisep_lowres2!=psisep_lowres: sys.exit('Inconsistent psisep_lowres: %7.4g, %7.4g' %(psisep_lowres,psisep_lowres2))

print "rmag_lowres", rmag_lowres
print "zmag_lowres", zmag_lowres
print "psiax_lowres", psiax_lowres
print "psisep_lowres", psisep_lowres
print "Bctr_lowres", Bctr_lowres
Rgrid_lowres = np.arange(nw_lowres)/float(nw_lowres-1)*rdim_lowres+rmin_lowres
print "rdim_lowres",rdim_lowres
print "rmin_lowres",rmin_lowres
#print "Rgrid_lowres",Rgrid_lowres
Zgrid_lowres = np.arange(nh_lowres)/float(nh_lowres-1)*zdim_lowres+(zmid_lowres-zdim_lowres/2.0)
print "zdim_lowres",zdim_lowres
print "zmid_lowres",zmid_lowres
#print "Zgrid_lowres",Zgrid_lowres

F_lowres=empty(nw_lowres,dtype=float)
p_lowres=empty(nw_lowres,dtype=float)
ffprime_lowres=empty(nw_lowres,dtype=float)
pprime_lowres=empty(nw_lowres,dtype=float)
qpsi_lowres=empty(nw_lowres,dtype=float)
psirz_1d_lowres=empty(nw_lowres*nh_lowres,dtype=float)
start_line=5
lines=range(nw_lowres/5)
if nw_lowres%5!=0: lines=range(nw_lowres/5+1)
for i in lines:
    n_entries=len(eqdsk_lowres[i+start_line])/entrylength
    F_lowres[i*5:i*5+n_entries]=[float(eqdsk_lowres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk_lowres[i+start_line])/entrylength
    p_lowres[i*5:i*5+n_entries]=[float(eqdsk_lowres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk_lowres[i+start_line])/entrylength
    ffprime_lowres[i*5:i*5+n_entries]=[float(eqdsk_lowres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk_lowres[i+start_line])/entrylength
    pprime_lowres[i*5:i*5+n_entries]=[float(eqdsk_lowres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

lines_twod=range(nw_lowres*nh_lowres/5)
if nw_lowres*nh_lowres%5!=0: lines_twod=range(nw_lowres*nh_lowres/5+1)
for i in lines_twod:
    n_entries=len(eqdsk_lowres[i+start_line])/entrylength
    psirz_1d_lowres[i*5:i*5+n_entries]=[float(eqdsk_lowres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1
psirz_lowres=psirz_1d_lowres.reshape(nh_lowres,nw_lowres)

for i in lines:
    n_entries=len(eqdsk_lowres[i+start_line])/entrylength
    qpsi_lowres[i*5:i*5+n_entries]=[float(eqdsk_lowres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

#linear grid of psi, on which all 1D fields are defined
linpsi_lowres=linspace(psiax_lowres,psisep_lowres,nw_lowres)
#create rho_tor_lowres grid
x_fine_lowres=linspace(psiax_lowres,psisep_lowres,nw_lowres*10)
phi_fine_lowres=empty((nw_lowres*10),dtype=float)
phi_fine_lowres[0]=0.
#spline of q for rhotor grid
interpol_order=3
q_spl_psi_lowres=US(linpsi_lowres,qpsi_lowres,k=interpol_order,s=1e-5)

for i in range(1,nw_lowres*10):
    x=x_fine_lowres[:i+1]
    y=q_spl_psi_lowres(x)
    phi_fine_lowres[i]=trapz(y,x)
rho_tor_fine_lowres=sqrt(phi_fine_lowres/phi_fine_lowres[-1])
rho_tor_spl_lowres=US(x_fine_lowres,rho_tor_fine_lowres,k=interpol_order,s=1e-5)
rho_tor_lowres=empty(nw_lowres,dtype=float)
for i in range(nw_lowres):
    rho_tor_lowres[i]=rho_tor_spl_lowres(linpsi_lowres[i])

rho_pol_fine_lowres = np.zeros(len(x_fine_lowres))
for i in range(len(x_fine_lowres)):
    rho_pol_fine_lowres[i] = sqrt((x_fine_lowres[i]-psiax_lowres)/(psisep_lowres-psiax_lowres))

if write_rhotor_rhopol_conversion:
    rt_rp_filename_lowres='rt_rp_%s' %filename_lowres
    rt_rp_file_lowres=open(rt_rp_filename_lowres,'w')
    rt_rp_file_lowres.write('# rho_tor_lowres          rho_pol_lowres\n')
    for i in range(len(x_fine_lowres)):
        rho_pol_lowres=sqrt((x_fine_lowres[i]-psiax_lowres)/(psisep_lowres-psiax_lowres))
        rt_rp_file_lowres.write('%16.8e %16.8e\n' %(rho_tor_fine_lowres[i],rho_pol_lowres))
    rt_rp_file_lowres.close()
    exit('\nWrote rhotor/rhopol conversion to %s.' %rt_rp_filename_lowres)

Z0_ind_lowres = np.argmin(np.abs(Zgrid_lowres-zmid_lowres))
psi_lowres=np.arange(nw_lowres)/float(nw_lowres-1)*(psisep_lowres-psiax_lowres)
#print "psi_lowres",psi_lowres
psi_midplane_lowres = psirz_lowres[Z0_ind_lowres,:]

#plt.plot(psi_midplane_lowres)
#plt.title(r'$\Psi$'+' midplane')
#plt.show()

#plt.plot(psi_midplane_lowres-psiax_lowres)
#plt.title(r'$\Psi$'+' midplane(adjusted)')
#plt.show()

#if show_plots:
    #plt.plot(psi_lowres/(psisep_lowres-psiax_lowres),qpsi_lowres,'x')
    #plt.xlabel(r'$\Psi/\Psi_0$')
    #plt.title(r'$q$')
    #plt.show()

    #plt.plot(psi_lowres/(psisep_lowres-psiax_lowres),p_lowres)
    #plt.xlabel(r'$\Psi/\Psi_0$')
    #plt.title(r'$P(Nm^{-2})$')
    #plt.show()

#plt.plot(psi_lowres,F_lowres)
#plt.xlabel(r'$\Psi$')
#plt.title(r'$RB_\phi(mT)$')
#plt.show()

    #plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),150)
    #plt.colorbar()
    #plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),levels = [1.0],color = 'black',linewidth = 2)
    #plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),levels = [0.9],color = 'black',linewidth = 2)
    #plt.xlabel('R(m)')
    #plt.ylabel('Z(m)')
    #plt.show()

#plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),50)
#plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),levels = [1.0],color = 'black', linewidth = 2)
    #plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),levels = [0.9,0.95,1.0],color = 'black', linewidth = 4)
    #plt.title('0.9,0.95,1.0')
    #plt.colorbar()
    #Zplot_lowres = np.zeros(nw_lowres)
    #for j in range(len(Zgrid_lowres)):
    #    Zplot_lowres[:] = Zgrid_lowres[j]
    #    plt.plot(Rgrid_lowres,Zplot_lowres,'.',color = 'black', markersize = 1)
    #plt.xlabel('R(m)')
    #plt.ylabel('Z(m)')
    #plt.show()

if write_binfo:
    #Calculate midplane B_pol
    B_pol_lowres = fd_d1_o4(psi_midplane_lowres,Rgrid_lowres)/Rgrid_lowres
    psi_norm_out_lowres = (psi_midplane_lowres-psiax_lowres)/(psisep_lowres-psiax_lowres)
    F_out_lowres = interp(psi_lowres/(psisep_lowres-psiax_lowres),F_lowres,psi_norm_out_lowres)
    q_out_lowres = interp(psi_lowres/(psisep_lowres-psiax_lowres),qpsi_lowres,psi_norm_out_lowres)
    rho_tor_out_lowres = interp(rho_pol_fine_lowres, rho_tor_fine_lowres, np.sqrt(psi_norm_out_lowres))
    f_lowres=open('Binfo_'+filename_lowres,'w')
    f_lowres.write('# Outer Midplane')
    f_lowres.write('# 1.R(m) 2.psi_norm 3.B_pol_lowres(T) 4.B_tor_lowres(T) 5.q 6.rho_tor_lowres\n')
    f_lowres.write('# R at magnetic axis = '+str(rmag_lowres)+'\n')
    f_lowres.write('# psisep_lowres - psiax_lowres = '+str(psisep_lowres-psiax_lowres)+'\n')
    Rmag_ind_lowres = np.argmin(abs(Rgrid_lowres - rmag_lowres))
    print "rmag_lowres",rmag_lowres
    print "Rmag_ind_lowres",Rmag_ind_lowres
    print "Rgrid_lowres[Rmag_ind_lowres]",Rgrid_lowres[Rmag_ind_lowres]
    temp = psi_norm_out_lowres
    temp[0:Rmag_ind_lowres] = 0
    psi_ind_sep_lowres = np.argmin(abs(temp-1.05))
    print "psi_ind_sep_lowres",psi_ind_sep_lowres
    B_tor_lowres = F_out_lowres / Rgrid_lowres
    np.savetxt(f_lowres,np.column_stack((Rgrid_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],psi_norm_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_pol_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_tor_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],q_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],rho_tor_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres])))
    f_lowres.close()
    #plt.plot(psi_norm_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],q_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres])
    #plt.plot(psi_lowres/(psisep_lowres-psiax_lowres),qpsi_lowres)
    #plt.plot(psi_norm_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],(Rgrid_lowres[Rmag_ind_lowres:psi_ind_sep_lowres]-rmag_lowres)*B_tor_lowres[Rmag_ind_lowres:psi_ind_sep_lowres]/(Rgrid_lowres[Rmag_ind_lowres:psi_ind_sep_lowres]*B_pol_lowres[Rmag_ind_lowres:psi_ind_sep_lowres]))
    #plt.show()
    
    #if show_plots:
        #plt.plot(Rgrid_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_pol_lowres[Rmag_ind_lowres:psi_ind_sep_lowres])
        #plt.title(r'$B_\theta$')
        #plt.xlabel('R(m)')
        #plt.show()

        #plt.plot(psi_norm_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_pol_lowres[Rmag_ind_lowres:psi_ind_sep_lowres]) 
        #plt.title(r'$B_\theta$')
        #plt.xlabel(r'$\Psi$')
        #plt.show()

eqdsk_highres=file_highres.readlines()
print 'Header: %s' %eqdsk_highres[0]
#set resolutions
nw_highres=int(eqdsk_highres[0].split()[-2]);nh_highres=int(eqdsk_highres[0].split()[-1])
pw_highres=(nw_highres/8/2)*2 #psi-width, number of flux surfaces around position of interest
print 'Resolution: %d x %d' %(nw_highres,nh_highres)

entrylength=16
#note: here rmin_highres is rleft from EFIT
try:
    rdim_highres,zdim_highres,rctr_highres,rmin_highres,zmid_highres=[float(eqdsk_highres[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_highres[1])/entrylength)]
except:
    entrylength=15
    try:
        rdim_highres,zdim_highres,rctr_highres,rmin_highres,zmid_highres=[float(eqdsk_highres[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_highres[1])/entrylength)]
    except:
        exit('Error reading EQDSK file, please check format!')

rmag_highres,zmag_highres,psiax_highres,psisep_highres,Bctr_highres=[float(eqdsk_highres[2][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_highres[2])/entrylength)]
dum_highres,psiax_highres2,dum_highres,rmag_highres2,dum_highres=[float(eqdsk_highres[3][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_highres[3])/entrylength)]
zmag_highres2,dum_highres,psisep_highres2,dum_highres,dum_highres=[float(eqdsk_highres[4][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk_highres[4])/entrylength)]
if rmag_highres!=rmag_highres2: sys.exit('Inconsistent rmag_highres: %7.4g, %7.4g' %(rmag_highres,rmag_highres2))
if psiax_highres2!=psiax_highres: sys.exit('Inconsistent psiax_highres: %7.4g, %7.4g' %(psiax_highres,psiax_highres2))
if zmag_highres!=zmag_highres2: sys.exit('Inconsistent zmag_highres: %7.4g, %7.4g' %(zmag_highres,zmag_highres2) )
if psisep_highres2!=psisep_highres: sys.exit('Inconsistent psisep_highres: %7.4g, %7.4g' %(psisep_highres,psisep_highres2))

print "rmag_highres", rmag_highres
print "zmag_highres", zmag_highres
print "psiax_highres", psiax_highres
print "psisep_highres", psisep_highres
print "Bctr_highres", Bctr_highres
Rgrid_highres = np.arange(nw_highres)/float(nw_highres-1)*rdim_highres+rmin_highres
print "rdim_highres",rdim_highres
print "rmin_highres",rmin_highres
#print "Rgrid_highres",Rgrid_highres
Zgrid_highres = np.arange(nh_highres)/float(nh_highres-1)*zdim_highres+(zmid_highres-zdim_highres/2.0)
print "zdim_highres",zdim_highres
print "zmid_highres",zmid_highres
#print "Zgrid_highres",Zgrid_highres

F_highres=empty(nw_highres,dtype=float)
p_highres=empty(nw_highres,dtype=float)
ffprime_highres=empty(nw_highres,dtype=float)
pprime_highres=empty(nw_highres,dtype=float)
qpsi_highres=empty(nw_highres,dtype=float)
psirz_1d_highres=empty(nw_highres*nh_highres,dtype=float)
start_line=5
lines=range(nw_highres/5)
if nw_highres%5!=0: lines=range(nw_highres/5+1)
for i in lines:
    n_entries=len(eqdsk_highres[i+start_line])/entrylength
    F_highres[i*5:i*5+n_entries]=[float(eqdsk_highres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk_highres[i+start_line])/entrylength
    p_highres[i*5:i*5+n_entries]=[float(eqdsk_highres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk_highres[i+start_line])/entrylength
    ffprime_highres[i*5:i*5+n_entries]=[float(eqdsk_highres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk_highres[i+start_line])/entrylength
    pprime_highres[i*5:i*5+n_entries]=[float(eqdsk_highres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

lines_twod=range(nw_highres*nh_highres/5)
if nw_highres*nh_highres%5!=0: lines_twod=range(nw_highres*nh_highres/5+1)
for i in lines_twod:
    n_entries=len(eqdsk_highres[i+start_line])/entrylength
    psirz_1d_highres[i*5:i*5+n_entries]=[float(eqdsk_highres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1
psirz_highres=psirz_1d_highres.reshape(nh_highres,nw_highres)

for i in lines:
    n_entries=len(eqdsk_highres[i+start_line])/entrylength
    qpsi_highres[i*5:i*5+n_entries]=[float(eqdsk_highres[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

#linear grid of psi, on which all 1D fields are defined
linpsi_highres=linspace(psiax_highres,psisep_highres,nw_highres)
#create rho_tor_highres grid
x_fine_highres=linspace(psiax_highres,psisep_highres,nw_highres*10)
phi_fine_highres=empty((nw_highres*10),dtype=float)
phi_fine_highres[0]=0.
#spline of q for rhotor grid
interpol_order=3
q_spl_psi_highres=US(linpsi_highres,qpsi_highres,k=interpol_order,s=1e-5)

for i in range(1,nw_highres*10):
    x=x_fine_highres[:i+1]
    y=q_spl_psi_highres(x)
    phi_fine_highres[i]=trapz(y,x)
rho_tor_fine_highres=sqrt(phi_fine_highres/phi_fine_highres[-1])
rho_tor_spl_highres=US(x_fine_highres,rho_tor_fine_highres,k=interpol_order,s=1e-5)
rho_tor_highres=empty(nw_highres,dtype=float)
for i in range(nw_highres):
    rho_tor_highres[i]=rho_tor_spl_highres(linpsi_highres[i])

rho_pol_fine_highres = np.zeros(len(x_fine_highres))
for i in range(len(x_fine_highres)):
    rho_pol_fine_highres[i] = sqrt((x_fine_highres[i]-psiax_highres)/(psisep_highres-psiax_highres))

if write_rhotor_rhopol_conversion:
    rt_rp_filename_highres='rt_rp_%s' %filename_highres
    rt_rp_file_highres=open(rt_rp_filename_highres,'w')
    rt_rp_file_highres.write('# rho_tor_highres          rho_pol_highres\n')
    for i in range(len(x_fine_highres)):
        rho_pol_highres=sqrt((x_fine_highres[i]-psiax_highres)/(psisep_highres-psiax_highres))
        rt_rp_file_highres.write('%16.8e %16.8e\n' %(rho_tor_fine_highres[i],rho_pol_highres))
    rt_rp_file_highres.close()
    exit('\nWrote rhotor/rhopol conversion to %s.' %rt_rp_filename_highres)

Z0_ind_highres = np.argmin(np.abs(Zgrid_highres-zmid_highres))
psi_highres=np.arange(nw_highres)/float(nw_highres-1)*(psisep_highres-psiax_highres)
#print "psi_highres",psi_highres
psi_midplane_highres = psirz_highres[Z0_ind_highres,:]

#plt.plot(psi_midplane_highres)
#plt.title(r'$\Psi$'+' midplane')
#plt.show()

#plt.plot(psi_midplane_highres-psiax_highres)
#plt.title(r'$\Psi$'+' midplane(adjusted)')
#plt.show()

if show_plots:
    plt.plot(psi_lowres/(psisep_lowres-psiax_lowres),qpsi_lowres, 'bo')
    plt.plot(psi_highres/(psisep_highres-psiax_highres),qpsi_highres,'r.')
    plt.xlabel('psi')
    plt.title('q')
    #plt.xlim((0.95,1.0))
    plt.show()

    plt.plot(psi_lowres/(psisep_lowres-psiax_lowres),p_lowres, 'bo')
    plt.plot(psi_highres/(psisep_highres-psiax_highres),p_highres, 'r.')
    plt.xlabel('psi_n')
    plt.title('P(N/m^2)')
    #plt.xlim((0.95,1.0))
    plt.show()

#plt.plot(psi_highres,F_highres)
#plt.xlabel(r'$\Psi$')
#plt.title(r'$RB_\phi(mT)$')
#plt.show()

if show_contour:
    #plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),150)
    #plt.colorbar()
    #plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),levels = [1.0],color = 'black',linewidth = 2)
    #plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),levels = [0.9],color = 'black',linewidth = 2)
    #plt.xlabel('R(m)')
    #plt.ylabel('Z(m)')
    #plt.show()

                 
#plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),50)
#plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),levels = [1.0],color = 'black', linewidth = 2)
    
    #plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),levels = [0.1,0.2,0.3],color = 'black', linewidth = 4)
    plt.contour(Rgrid_lowres,Zgrid_lowres,(psirz_lowres[:]-psiax_lowres)/(psisep_lowres-psiax_lowres),levels = [0.9, 0.95, 0.999],color = 'black', linewidth = 4)
    plt.title('low res')
    #plt.colorbar()
    Zplot_lowres = np.zeros(nw_lowres)
    for j in range(len(Zgrid_lowres)):
        Zplot_lowres[:] = Zgrid_lowres[j]
        plt.plot(Rgrid_lowres,Zplot_lowres,'.',color = 'black', markersize = 1)
    	plt.xlabel('R(m)')
    	plt.ylabel('Z(m)')
    #plt.show()

#    plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),levels = [0.1,0.2,0.3],color = 'black', linewidth = 4)
    plt.contour(Rgrid_highres,Zgrid_highres,(psirz_highres[:]-psiax_highres)/(psisep_highres-psiax_highres),levels = [0.9, 0.95, 1.0],color = 'black', linewidth = 4)
    plt.title('high res')
    plt.colorbar()
    #Zplot_highres = np.zeros(nw_highres)
    #for j in range(len(Zgrid_highres)):
    #    Zplot_highres[:] = Zgrid_highres[j]
    #    plt.plot(Rgrid_highres,Zplot_highres,'.',color = 'black', markersize = 1)
    #plt.xlabel('R(m)')
    #plt.ylabel('Z(m)')
    plt.show()

if write_binfo:
    #Calculate midplane B_pol
    B_pol_highres = fd_d1_o4(psi_midplane_highres,Rgrid_highres)/Rgrid_highres
    psi_norm_out_highres = (psi_midplane_highres-psiax_highres)/(psisep_highres-psiax_highres)
    F_out_highres = interp(psi_highres/(psisep_highres-psiax_highres),F_highres,psi_norm_out_highres)
    q_out_highres = interp(psi_highres/(psisep_highres-psiax_highres),qpsi_highres,psi_norm_out_highres)
    rho_tor_out_highres = interp(rho_pol_fine_highres, rho_tor_fine_highres, np.sqrt(psi_norm_out_highres))
    f_highres=open('Binfo_'+filename_highres,'w')
    f_highres.write('# Outer Midplane')
    f_highres.write('# 1.R(m) 2.psi_norm 3.B_pol_highres(T) 4.B_tor_highres(T) 5.q 6.rho_tor_highres\n')
    f_highres.write('# R at magnetic axis = '+str(rmag_highres)+'\n')
    f_highres.write('# psisep_highres - psiax_highres = '+str(psisep_highres-psiax_highres)+'\n')
    Rmag_ind_highres = np.argmin(abs(Rgrid_highres - rmag_highres))
    print "rmag_highres",rmag_highres
    print "Rmag_ind_highres",Rmag_ind_highres
    print "Rgrid_highres[Rmag_ind_highres]",Rgrid_highres[Rmag_ind_highres]
    temp = psi_norm_out_highres
    temp[0:Rmag_ind_highres] = 0
    psi_ind_sep_highres = np.argmin(abs(temp-1.05))
    print "psi_ind_sep_highres",psi_ind_sep_highres
    B_tor_highres = F_out_highres / Rgrid_highres
    np.savetxt(f_highres,np.column_stack((Rgrid_highres[Rmag_ind_highres:psi_ind_sep_highres],psi_norm_out_highres[Rmag_ind_highres:psi_ind_sep_highres],B_pol_highres[Rmag_ind_highres:psi_ind_sep_highres],B_tor_highres[Rmag_ind_highres:psi_ind_sep_highres],q_out_highres[Rmag_ind_highres:psi_ind_sep_highres],rho_tor_out_highres[Rmag_ind_highres:psi_ind_sep_highres])))
    f_highres.close()
    #plt.plot(psi_norm_out_highres[Rmag_ind_highres:psi_ind_sep_highres],q_out_highres[Rmag_ind_highres:psi_ind_sep_highres])
    #plt.plot(psi_highres/(psisep_highres-psiax_highres),qpsi_highres)
    #plt.plot(psi_norm_out_highres[Rmag_ind_highres:psi_ind_sep_highres],(Rgrid_highres[Rmag_ind_highres:psi_ind_sep_highres]-rmag_highres)*B_tor_highres[Rmag_ind_highres:psi_ind_sep_highres]/(Rgrid_highres[Rmag_ind_highres:psi_ind_sep_highres]*B_pol_highres[Rmag_ind_highres:psi_ind_sep_highres]))
    #plt.show()

    if show_plots:
	#plt.plot(Rgrid_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_pol_lowres[Rmag_ind_lowres:psi_ind_sep_lowres])
        #plt.plot(Rgrid_highres[Rmag_ind_highres:psi_ind_sep_highres],B_pol_highres[Rmag_ind_highres:psi_ind_sep_highres])
        #plt.title(r'$B_\theta$')
	#plt.xlabel('R(m)')
        #plt.show()

	plt.plot(psi_norm_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_pol_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],'bx')
        plt.plot(psi_norm_out_highres[Rmag_ind_highres:psi_ind_sep_highres],B_pol_highres[Rmag_ind_highres:psi_ind_sep_highres],'r.')
        plt.title('B_pol')
        plt.xlabel('psi_n')
	#plt.xlim((0.95,1.0))
        plt.show()

	plt.plot(psi_norm_out_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],B_tor_lowres[Rmag_ind_lowres:psi_ind_sep_lowres],'bx')
        plt.plot(psi_norm_out_highres[Rmag_ind_highres:psi_ind_sep_highres],B_tor_highres[Rmag_ind_highres:psi_ind_sep_highres],'r.')
        plt.title('B_tor')
        plt.xlabel('psi_n')
	#plt.xlim((0.95,1.0))
        plt.show()


             
