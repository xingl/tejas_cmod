#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import math
import optparse as op

parser = op.OptionParser()
(options,args) = parser.parse_args()


file_name_1 = args[0]
file_name_2 = args[1]



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
                

f_1=open(file_name_1,'r')
data = f_1.read()
data_lines = data.split('\n')

rbs_string = []
for line in data_lines:
    if 'RBC' in line or 'RBS' in line or 'ZBC' in line or 'ZBS' in line:
        rbs_string.append(line)

rbs_1 = np.empty(0)
rbc_1 = np.empty(0)
zbs_1 = np.empty(0)
zbc_1 = np.empty(0)

for i in rbs_string:
    num = extract_coeff_num(i)
    this_rbs = extract_starbs(i,'rbs')
    if this_rbs == this_rbs:
        rbs_1 = np.append(rbs_1,this_rbs)
    this_rbc = extract_starbs(i,'rbc')
    if this_rbc == this_rbc:
        rbc_1 = np.append(rbc_1,this_rbc)
    this_zbs = extract_starbs(i,'zbs')
    if this_zbs == this_zbs:
        zbs_1 = np.append(zbs_1,this_zbs)
    this_zbc = extract_starbs(i,'zbc')
    if this_zbc == this_zbc:
        zbc_1 = np.append(zbc_1,this_zbc)

num_1 = np.arange(len(rbs_1),dtype='int')

f_2=open(file_name_2,'r')
data = f_2.read()
data_lines = data.split('\n')

rbs_string = []
for line in data_lines:
    if 'RBC' in line or 'RBS' in line or 'ZBC' in line or 'ZBS' in line:
        rbs_string.append(line)

rbs_2 = np.empty(0)
rbc_2 = np.empty(0)
zbs_2 = np.empty(0)
zbc_2 = np.empty(0)

for i in rbs_string:
    num = extract_coeff_num(i)
    this_rbs = extract_starbs(i,'rbs')
    if this_rbs == this_rbs:
        rbs_2 = np.append(rbs_2,this_rbs)
    this_rbc = extract_starbs(i,'rbc')
    if this_rbc == this_rbc:
        rbc_2 = np.append(rbc_2,this_rbc)
    this_zbs = extract_starbs(i,'zbs')
    if this_zbs == this_zbs:
        zbs_2 = np.append(zbs_2,this_zbs)
    this_zbc = extract_starbs(i,'zbc')
    if this_zbc == this_zbc:
        zbc_2 = np.append(zbc_2,this_zbc)

num_2 = np.arange(len(rbs_2),dtype='int')

plt.plot(num_1,rbs_1,label='rbs_1')
plt.plot(num_1,rbc_1,label='rbc_1')
plt.plot(num_1,zbs_1,label='zbs_1')
plt.plot(num_1,zbc_1,label='zbc_1')
plt.plot(num_2,rbs_2,label='rbs_2')
plt.plot(num_2,rbc_2,label='rbc_2')
plt.plot(num_2,zbs_2,label='zbs_2')
plt.plot(num_2,zbc_2,label='zbc_2')
plt.legend()
plt.xlabel('i')
plt.show()

r_1 = np.empty(0)
z_1 = np.empty(0)
theta_1 = np.empty(0)

for i in range(10000):
    this_theta = i*2*math.pi/10000
    this_r = 0.0
    this_z = 0.0
    for j in range(len(num_1)):
        arg = j*this_theta
        this_r = this_r + rbc_1[j]*math.cos(arg) + rbs_1[j]*math.sin(arg)
        this_z = this_z + zbc_1[j]*math.cos(arg) + zbs_1[j]*math.sin(arg)
    r_1 = np.append(r_1, this_r)
    z_1 = np.append(z_1, this_z)
    theta_1 = np.append(theta_1, this_theta)

r_2 = np.empty(0)
z_2 = np.empty(0)
theta_2 = np.empty(0)

for i in range(10000):
    this_theta = i*2*math.pi/10000
    this_r = 0.0
    this_z = 0.0
    for j in range(len(num_2)):
        arg = j*this_theta
        this_r = this_r + rbc_2[j]*math.cos(arg) + rbs_2[j]*math.sin(arg)
        this_z = this_z + zbc_2[j]*math.cos(arg) + zbs_2[j]*math.sin(arg)
    r_2 = np.append(r_2, this_r)
    z_2 = np.append(z_2, this_z)
    theta_2 = np.append(theta_2, this_theta)

plt.plot(r_1,z_1)
plt.plot(r_2,z_2)
plt.xlabel('r')
plt.ylabel('z')
plt.show()

