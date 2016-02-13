import numpy as np
from interp import *
import matplotlib.pyplot as plt


def read_cmod_pfile(p_file_name,shift_Te=False,shift_Ti=True,add_impurity=True):

    impurity_charge = 5.0

    f=open(p_file_name,'r')
    data = f.read()
    f.close()

    sdata = data.split('\n')
    nr = int(sdata[0].split()[0])
    print "p-file resolution: nr = ", nr

    # ne is electron density profile, dne is gradient of ne
    ne = np.empty(0)
    dne = np.empty(0)
    ni = np.empty(0)
    dni = np.empty(0)
    te = np.empty(0)
    dte = np.empty(0)
    ti = np.empty(0)
    dti = np.empty(0)

    # psipne is the grid of psi_pol on which ne&dne above is recorded
    psipne = np.empty(0)
    psipni = np.empty(0)
    psipte = np.empty(0)
    psipti = np.empty(0)

    for i in np.array(range(nr)):
        temp = sdata[i+1].split()
        psipne = np.append(psipne,float(temp[0]))
        ne = np.append(ne,float(temp[1]))
        dne = np.append(dne,float(temp[2]))
        temp = sdata[nr+i+2].split()
        psipte = np.append(psipte,float(temp[0]))
        te = np.append(te,float(temp[1]))
        dte = np.append(dte,float(temp[2]))
        temp = sdata[2*nr+i+3].split()
        psipni = np.append(psipni,float(temp[0]))
        ni = np.append(ni,float(temp[1]))
        dni = np.append(dni,float(temp[2]))
        temp = sdata[3*nr+i+4].split()
        psipti = np.append(psipti,float(temp[0]))
        ti = np.append(ti,float(temp[1]))
        dti = np.append(dti,float(temp[2]))

    if shift_Te:
       psipte = psipte + 0.005
    if shift_Ti:
       psipti = psipti + 0.005

    #plt.plot(psipne,'x')
    #plt.plot(psipni,'r.')
    #plt.show()

    #psi0 is normalized psi_pol in [0,1] with 1000 points
    #quantities with _out are interpolated on psi0 grid 
    #print "length of arrays",len(ne),len(ni),len(te),len(ti)
    psi0 = np.arange(1000)/999.0
    #psi0 = np.linspace(0.0,1.0,len(ne))
    ne_out = interp(psipne,ne,psi0)
    te_out = interp(psipte,te,psi0)
    ni_out = interp(psipni,ni,psi0)
    ti_out = interp(psipti,ti,psi0)
    dne_out = interp(psipne,dne,psi0)
    dte_out = interp(psipte,dte,psi0)
    dni_out = interp(psipni,dni,psi0)
    dti_out = interp(psipti,dti,psi0)
    
    if add_impurity:
       nz_out = np.empty(len(ne_out))
       for i in range(len(ne_out)):
           nz_out[i] = (ne_out[i]-ni_out[i])/impurity_charge


    # output quantities are psi_pol, ne, te, ni, ti
    # grid: even psi_pol_norm
    # resolution: 1000
    return psi0, ne_out, te_out, ni_out, ti_out, nz_out


    
def read_cmod_pfile_raw(p_file_name):

    impurity_charge = 5.0

    f=open(p_file_name,'r')
    data = f.read()
    f.close()

    sdata = data.split('\n')
    nr = int(sdata[0].split()[0])
    print "p file resolution: nr = ", nr

    # ne is electron density profile, dne is gradient of ne
    ne = np.empty(0)
    dne = np.empty(0)
    ni = np.empty(0)
    dni = np.empty(0)
    te = np.empty(0)
    dte = np.empty(0)
    ti = np.empty(0)
    dti = np.empty(0)
    # psipne is the grid of psi_pol on which ne&dne above is recorded
    psipne = np.empty(0)
    psipni = np.empty(0)
    psipte = np.empty(0)
    psipti = np.empty(0)

    for i in np.array(range(nr)):
        temp = sdata[i+1].split()
        psipne = np.append(psipne,float(temp[0]))
        ne = np.append(ne,float(temp[1]))
        dne = np.append(dne,float(temp[2]))
        temp = sdata[nr+i+2].split()
        psipte = np.append(psipte,float(temp[0]))
        te = np.append(te,float(temp[1]))
        dte = np.append(dte,float(temp[2]))
        temp = sdata[2*nr+i+3].split()
        psipni = np.append(psipni,float(temp[0]))
        ni = np.append(ni,float(temp[1]))
        dni = np.append(dni,float(temp[2]))
        temp = sdata[3*nr+i+4].split()
        psipti = np.append(psipti,float(temp[0]))
        ti = np.append(ti,float(temp[1]))
        dti = np.append(dti,float(temp[2]))


    return psipne,ne,psipte,te,psipni,ni,psipti,ti
