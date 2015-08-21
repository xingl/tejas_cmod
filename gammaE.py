from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
import re

parser=op.OptionParser()
options,args=parser.parse_args()

efit_filename=args[0]
efit_file=open(efit_filename,'r')

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

eqdsk=efit_file.readlines()
nw=int(eqdsk[0].split()[-2]);nh=int(eqdsk[0].split()[-1])
pw=(nw/8/2)*2 #psi-width, number of flux surfaces around position of interest
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

Rgrid = np.arange(nw)/float(nw-1)*rdim+rmin
Zgrid = np.arange(nh)/float(nh-1)*zdim+(zmid-zdim/2.0)

F=empty(nw,dtype=float)
p=empty(nw,dtype=float)
ffprime=empty(nw,dtype=float)
pprime=empty(nw,dtype=float)
psirz_1d=empty(nw*nh,dtype=float)
qpsi=empty(nw,dtype=float)
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



#linear grid of psi, on which all 1D fields are defined
linpsi=linspace(psiax,psisep,nw)
#create rho_tor grid
x_fine=linspace(psiax,psisep,nw*10)
#print "x_fine", x_fine[0:10]
phi_fine=empty((nw*10),dtype=float)
phi_fine[0]=0.
interpol_order=3
q_spl_psi=US(linpsi,qpsi,k=interpol_order,s=1e-5)

for i in range(1,nw*10):
    x=x_fine[:i+1]
    y=q_spl_psi(x)
    phi_fine[i]=trapz(y,x)

#plt.plot(x_fine, phi_fine, 'rx')
#plt.show()

rho_tor_fine=sqrt(phi_fine/phi_fine[-1])
#plt.plot(x_fine, rho_tor_fine, 'rx')
#plt.show()
q_fine = q_spl_psi(x_fine)

rho_tor_spl=US(x_fine,rho_tor_fine,k=interpol_order,s=1e-5)
rho_tor=empty(nw,dtype=float)

for i in range(nw):
    rho_tor[i]=rho_tor_spl(linpsi[i])

rho_pol_fine = np.zeros(len(x_fine))
for i in range(len(x_fine)):
    rho_pol_fine[i] = sqrt((x_fine[i]-psiax)/(psisep-psiax))


Z0_ind = np.argmin(np.abs(Zgrid-zmid))
psi=np.arange(nw)/float(nw-1)*(psisep-psiax)
#print "psi",psi
psi_midplane = psirz[Z0_ind,:]

Bp_evenR = fd_d1_o4(psi_midplane,Rgrid)/Rgrid
psi_norm_evenR = (psi_midplane-psiax)/(psisep-psiax)
F_evenR = interp(psi/(psisep-psiax),F,psi_norm_evenR)
p_evenR = interp(psi/(psisep-psiax),p,psi_norm_evenR)
rho_tor_evenR = interp(rho_pol_fine, rho_tor_fine, np.sqrt(psi_norm_evenR))

f = open('pb.dat','w')
f.write('# Outer Midplane')
f.write('# 1.R(m) 2.psi_norm 3.B_pol(T) 4.B_tor(T) 5.q 6.Pressure 7.rho_tor \n')
f.write('# R at magnetic axis = '+str(rmag)+'\n')
f.write('# psisep - psiax = '+str(psisep-psiax)+'\n')
Rmag_ind = np.argmin(abs(Rgrid - rmag))
print "rmag",rmag
print "Rmag_ind",Rmag_ind
print "Rgrid[Rmag_ind]",Rgrid[Rmag_ind]
print "psi_norm_evenR[Rmag_ind]",psi_norm_evenR[Rmag_ind]
print "psiax", psiax

temp = np.copy(psi_norm_evenR)
temp[0:Rmag_ind] = 0
psi_ind_sep = np.argmin(abs(temp-1.05))
psisep_ind = np.argmin(abs(temp-1.0)) 
print "psisep_ind",psisep_ind
print "psi_norm_evenR[psisep_ind]", psi_norm_evenR[psisep_ind]
print "psi_midplane[psisep_ind]", psi_midplane[psisep_ind]
print "psisep", psisep

Bt_evenR = F_evenR / Rgrid
p_evenR = p_evenR/np.max(p_evenR)
q_evenR = interp(psi/(psisep-psiax),qpsi,psi_norm_evenR)
np.savetxt(f,np.column_stack((Rgrid[Rmag_ind:psisep_ind],psi_norm_evenR[Rmag_ind:psisep_ind],Bp_evenR[Rmag_ind:psisep_ind],Bt_evenR[Rmag_ind:psisep_ind],q_evenR[Rmag_ind:psisep_ind],p_evenR[Rmag_ind:psisep_ind],rho_tor_evenR[Rmag_ind:psisep_ind])))
f.close()


n=3000
evenpsi = linspace(psi_norm_evenR[Rmag_ind],psi_norm_evenR[psisep_ind],n)

print "psi_norm_evenR[psisep_ind]", psi_norm_evenR[psisep_ind]
print "evenpsi[-1]", evenpsi[-1]

R_evenpsi = interp(psi_norm_evenR[Rmag_ind:psisep_ind],
	    Rgrid[Rmag_ind:psisep_ind],evenpsi) 
print "Rgrid[psisep_ind]", Rgrid[psisep_ind]
print "R_evenpsi[-1]", R_evenpsi[-1]

Bp_evenpsi = interp(Rgrid,Bp_evenR,R_evenpsi) 
F_evenpsi = interp(psi/(psisep-psiax),F,evenpsi) 
Bt_evenpsi = F_evenpsi / R_evenpsi 
p_evenpsi = interp(psi/(psisep-psiax),p,evenpsi) 
p_evenpsi = p_evenpsi/np.max(p_evenpsi) 
q_evenpsi = interp(psi/(psisep-psiax),qpsi,evenpsi)

rhot_evenpsi = interp(rho_pol_fine,rho_tor_fine,np.sqrt(evenpsi))

f2 = open('evenpsi.dat','w') 
f2.write('# 1.R_evenpsi 2.evenpsi 3.B_p 4.B_t 5.q_out 6.P 7.rhot_evenpsi\n') 
np.savetxt(f2,np.column_stack((R_evenpsi,evenpsi,Bp_evenpsi,Bt_evenpsi,q_evenpsi,p_evenpsi,rhot_evenpsi)))
f2.close()


plt.plot(Rgrid[Rmag_ind:psisep_ind],psi_norm_evenR[Rmag_ind:psisep_ind],'x')
plt.plot(R_evenpsi,evenpsi,'r.')
plt.xlabel('R')
plt.ylabel('psi_p')
plt.show()

plt.plot(psi_norm_evenR[Rmag_ind:psisep_ind],Bp_evenR[Rmag_ind:psisep_ind],'x')
plt.plot(evenpsi,Bp_evenpsi,'r.')
plt.xlabel('psi_p')
plt.ylabel('B_p')
plt.show()

plt.plot(psi_norm_evenR[Rmag_ind:psisep_ind],Bt_evenR[Rmag_ind:psisep_ind],'x')
plt.plot(evenpsi,Bt_evenpsi,'r.')
plt.xlabel('psi_p')
plt.ylabel('B_t')
plt.show()

plt.plot(psi_norm_evenR[Rmag_ind:psisep_ind],q_evenR[Rmag_ind:psisep_ind],'x')
plt.plot(evenpsi,q_evenpsi,'r.')
plt.xlabel('psi_p')
plt.ylabel('q')
plt.show()

plt.plot(psi_norm_evenR[Rmag_ind:psisep_ind],p_evenR[Rmag_ind:psisep_ind],'x')
plt.plot(evenpsi,p_evenpsi,'r.')
plt.xlabel('psi_p')
plt.ylabel('p')
plt.show()

plt.plot(psi_norm_evenR[Rmag_ind:psisep_ind],rho_tor_evenR[Rmag_ind:psisep_ind],'x')
plt.plot(evenpsi,rhot_evenpsi,'r.')
plt.xlabel('psi_p')
plt.ylabel('rho_t')
plt.show()


iterdb_filename=args[1]
iterdb_file=open(iterdb_filename,'r')

data_in=iterdb_file.read()
data_linesplit=data_in.split('\n')

keep_going=1
i=0
while keep_going:
    test=re.search(';-# OF X PTS',data_linesplit[i])
    if test:
        num=data_linesplit[i].split()[0]
        num=float(num)
        num=int(num)
        print "number of points:",num
        keep_going=(1==2)
    if i == len(data_linesplit):
        keep_going=(1==2)
    i=i+1

lnum=0
while len(data_linesplit)-lnum > 10:
    sec_num_lines = num/6
    if num % 6 != 0:
        sec_num_lines += 1
    keep_going=1
    while keep_going:
        test=re.search('-DEPENDENT VARIABLE LABEL',data_linesplit[lnum])
        #test2=re.search('INDEPENDENT',data_linesplit[lnum])
        if test :
            quantity=data_linesplit[lnum].split()[0]
            units=data_linesplit[lnum].split()[1]
        test2=re.search('DATA FOLLOW',data_linesplit[lnum])
        if test2:
            keep_going=(1==2)
        lnum=lnum+1

    if quantity=='NM1':
       print "Reading :",quantity
       rhot1=np.empty(0)
       lnum0 = lnum
       for j in range(lnum0,lnum0+sec_num_lines):
           for k in range(6):
               str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
               if(re.search('e',str_temp)):
                   temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                   rhot1=np.append(rhot1,temp)
           lnum=lnum+1
       
       lnum=lnum+1

       arr1=np.empty(0)
       lnum0 = lnum
       for j in range(lnum0,lnum0+sec_num_lines):
           for k in range(6):
               str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
               if(re.search('e',str_temp)):
                  temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                  arr1=np.append(arr1,temp)
           lnum=lnum+1

    if quantity=='TI':
       print "Reading :",quantity
       rhot2=np.empty(0)
       lnum0 = lnum
       for j in range(lnum0,lnum0+sec_num_lines):
           for k in range(6):
               str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
               if(re.search('e',str_temp)):
                   temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                   rhot2=np.append(rhot2,temp)
           lnum=lnum+1

       lnum=lnum+1

       arr2=np.empty(0)
       lnum0 = lnum
       for j in range(lnum0,lnum0+sec_num_lines):
           for k in range(6):
               str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
               if(re.search('e',str_temp)):
                  temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                  arr2=np.append(arr2,temp)
           lnum=lnum+1
    else:
       lnum = lnum + 2*sec_num_lines + 1
	
vout=np.empty((len(arr1),4),dtype='float')
vout[:,0]=rhot1
vout[:,1]=arr1/np.max(arr1)
vout[:,2]=rhot2
vout[:,3]=arr2/np.max(arr2)
f=open('nt.dat','w')
f.write('#'+iterdb_filename+'\n'+'#1.rhot 2.density 3.rhot 4.temp_i'+ '\n')
np.savetxt(f,vout)
f.close()

dens = vout[:,1]
tempi= vout[:,3]

plt.plot(rhot_evenpsi,p_evenpsi,'x')
plt.plot(rhot1,dens*tempi,'r.')
plt.xlabel('rho_t')
plt.ylabel('P')
plt.title('blue from EFIT file, red from ITERDB file')
plt.show()

#plt.plot(rho_tor_evenR[Rmag_ind:psisep_ind],p_evenR[Rmag_ind:psisep_ind],'b.')
#plt.plot(rhot1,dens*tempi,'rx')
#plt.show()

dens_evenpsi = interp(rhot1,dens,rhot_evenpsi)
tempi_evenpsi = interp(rhot1,tempi,rhot_evenpsi)

plt.plot(rhot1,dens,'x')
plt.plot(rhot_evenpsi,dens_evenpsi,'r.')
plt.xlabel('rho_t')
plt.ylabel('n')
plt.show()

plt.plot(rhot1,tempi,'x')
plt.plot(rhot_evenpsi,tempi_evenpsi,'r.')
plt.xlabel('rho_t')
plt.ylabel('Ti')
plt.show()


Er = fd_d1_o4(dens_evenpsi*tempi_evenpsi,evenpsi)/dens_evenpsi
B_tot = sqrt(Bp_evenpsi**2+Bt_evenpsi**2)
gammaE = fd_d1_o4(Er/B_tot,evenpsi)

plt.plot(rhot_evenpsi,gammaE)
plt.ylabel('gamma_E')
plt.xlabel('rho_t (at even psi_p)')
plt.show()

evenrhot = linspace(rhot_evenpsi[0],rhot_evenpsi[-1],n)
q_evenrhot = interp(rhot_evenpsi,q_evenpsi,evenrhot)

plt.plot(evenrhot,q_evenrhot,'x')
plt.plot(rhot_evenpsi,q_evenpsi,'r.')
plt.xlabel('rhot')
plt.ylabel('q')
plt.show()

shat = evenrhot*fd_d1_o4(q_evenrhot,evenrhot)/q_evenrhot
plt.plot(evenrhot,shat)
plt.xlabel('even rho_t')
plt.ylabel('shat')
plt.show()



