"""
Read the first column ( square root of psi_tor ) and 
the second column ( psi_pol ) from the 'jbsL_  /rbsProfs'.

Read the second column ( psi_pol ) and the nineth 
column ( ballooning instability indicator ) from the 
file 'dcon.out' in the directory 'jbsL_   /VM_new_1/'.

Use the function 'interp' to interpolate the relation 
between square root of psi_tor and psi_pol and translate 
the values of psi_pol into square root of psi_tor. 

Plot ballooning instability indicator as a function of 
square root of psi_tor.
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import optparse as op

parser = op.OptionParser()
(options,args) = parser.parse_args()

if len(args) != 1:
        exit("""Please include a VMEC input file name as an argument.\n """)

plot_title = args[0]

#file_name = 'rbsProfs'
#f = open(file_name,'r')
#data = f.read()
#f.close()
#data_lines = data.split('\n')

#tor_psi_sqrt = np.empty(0)
#pol_psi = np.empty(0)

#for line in data_lines:
#        if line != '':
#                if line[0] != '#':
#                        line_split = line.split()
#                        tor_psi_sqrt = np.append(tor_psi_sqrt, float(line_split[0]))
#                        pol_psi = np.append(pol_psi, float(line_split[1]))

#f=open('pol_psi_tor_psi_sqrt.temp','w')
#f.write('#'+file_name+'\n')
#f.write('# tor_psi_sqrt  pol_psi'+'\n')
#np.savetxt(f,np.column_stack((tor_psi_sqrt, pol_psi)))
#f.close()

data = np.genfromtxt('rbsProfs')
tor_psi_sqrt = data[:,0]
pol_psi = data[:,1]

#f=open('fromrbsProfs2.temp','w')
#f.write('#'+file_name+'\n')
#f.write('# tor_psi_sqrt  pol_psi'+'\n')
#np.savetxt(f,np.column_stack((tor_psi_sqrt2, pol_psi2)))
#f.close()

file_name = 'dcon.out'
f = open(file_name,'r')
data = f.read()
f.close()
data_lines = data.split('\n')

pol_psi_fac = np.empty(0)
balloon = np.empty(0)
n = np.empty(0)

def isfloat(s):
        try:
                float(s)
                return True
        except ValueError:
                return False

for line in data_lines:
                if 'ca1' in line:
                        if len(n) == 0:
                                n = np.append(n, data_lines.index(line))
                        if len(n) == 1:
                                n = np.append(n, data_lines.index(line, int(n[0])+1))

for i in range(int(n[0])+1,int(n[1])+1):
        if data_lines[i] != '':
                line_split = data_lines[i].split()
                if isfloat(line_split[8]):
                        pol_psi_fac = np.append(pol_psi_fac,float(line_split[1]))
                        balloon = np.append(balloon,float(line_split[8]))

#f=open('balloon_pol.temp','w')
#f.write('#'+file_name+'\n')
#f.write('# pol balloon'+'\n')
#np.savetxt(f,np.column_stack((pol_psi_fac,balloon)))
#f.close()


def interp(xin,yin,xnew):
	rho_tck = interpolate.splrep(xin,yin)
	yout = interpolate.splev(xnew,rho_tck,der=0)
	
	return yout

new_tor_psi_sqrt = interp(pol_psi,tor_psi_sqrt,pol_psi_fac)

plt.plot(new_tor_psi_sqrt,balloon)
plt.title(plot_title)
plt.xlabel("sqrt_psi_tor")
plt.show()
