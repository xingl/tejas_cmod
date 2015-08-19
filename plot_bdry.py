"""
Plot the shape of the last closed magnetic surface from 
VMEC input files.

Takes the full input file name as an argument.

Fourier Transform the coefficients in front of sin and 
cos of R (RBS, RBC) and Z (ZBS, ZBC) in to 
r(theta) and z(theta).

Plot r,z.
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import optparse as op

parser = op.OptionParser()
(options,args) = parser.parse_args()

if len(args) != 1:
	exit("""Please include a VMEC input file name as an argument.\n """)

file_name = args[0]

#file_name = 'input_n70_0.999.akima.1120907032.01012'
f=open(file_name,'r')
data = f.read()
data_lines = data.split('\n')

rbs_string = []
for line in data_lines:
    if 'RBC' in line or 'RBS' in line or 'ZBC' in line or 'ZBS' in line:
        rbs_string.append(line)

def extract_coeff_num(line):
    num_out = '' 
    keep_going = True
    line2 = line.replace(' ','')
    j = 1
    for i in range(len(line2)):
        if line2[i]==',':
            while keep_going:
                if line2[i+j] != ')':
                    num_out += line2[i+j]
                    j += 1
                else:
                    keep_going = False
                
    return int(float(num_out))

def extract_starbs(line,star_bs):
    if star_bs in line.lower():
        line_split = line.split()
        for i in range(len(line_split)):
            if star_bs in line_split[i].lower():
                rbs_out = float(line_split[i+2])
    else:
        rbs_out = np.nan
    return rbs_out
                

#print rbs_string
rbs = np.empty(0)
rbc = np.empty(0)
zbs = np.empty(0)
zbc = np.empty(0)
for i in rbs_string:
    num = extract_coeff_num(i)
    this_rbs = extract_starbs(i,'rbs')
    if this_rbs == this_rbs:
        rbs = np.append(rbs,this_rbs)
    this_rbc = extract_starbs(i,'rbc')
    if this_rbc == this_rbc:
        rbc = np.append(rbc,this_rbc)
    this_zbs = extract_starbs(i,'zbs')
    if this_zbs == this_zbs:
        zbs = np.append(zbs,this_zbs)
    this_zbc = extract_starbs(i,'zbc')
    if this_zbc == this_zbc:
        zbc = np.append(zbc,this_zbc)

num = np.arange(len(rbs),dtype='int')
#f=open('rbs_rbc_zbs_zbc.temp','w')
#f.write('#'+file_name+'\n')
#f.write('# num rbc rbs zbc zbs'+'\n')
#np.savetxt(f,np.column_stack((num,rbc,rbs,zbc,zbs)))
#f.close()

r = np.empty(0)
z = np.empty(0)
theta = np.empty(0)

for i in range(2000):
    this_theta = i*2*math.pi/360.0
    this_r = 0.0
    this_z = 0.0
    for j in range(len(num)):
        arg = j*this_theta
        this_r = this_r + rbc[j]*math.cos(arg) + rbs[j]*math.sin(arg)
        this_z = this_z + zbc[j]*math.cos(arg) + zbs[j]*math.sin(arg)
    r = np.append(r, this_r)
    z = np.append(z, this_z)
    theta = np.append(theta, this_theta)

#f = open('theta_r_z.temp','w')
#f.write('#'+file_name+'\n')
#f.write('# num theta r z '+'\n')
#np.savetxt(f,np.column_stack((num,theta,r,z)))
#f.close()

plt.plot(r,z)
plt.title(file_name)
plt.xlabel('r')
plt.ylabel('z')
plt.show()
#plt.savefig(file_name+'.pdf', format='pdf')
