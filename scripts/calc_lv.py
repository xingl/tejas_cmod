#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from interp import *

parser = op.OptionParser()
options,args = parser.parse_args()
c_buffer_x = args[0]
buffer_size = args[1]
l_buffer_x = float(c_buffer_x) - 0.35*float(buffer_size)
r_buffer_x = float(c_buffer_x) + 0.4*float(buffer_size)

pdata=np.genfromtxt('p_info.dat')
rhot=pdata[:,0]
te=pdata[:,2]

rhot_fine = linspace(rhot[0],rhot[-1],10*len(rhot))
te_fine = interp(rhot,te,rhot_fine)

l_ind = np.argmin(abs(rhot_fine - float(l_buffer_x)))
c_ind = np.argmin(abs(rhot_fine - float(c_buffer_x)))

te_l = te_fine[l_ind]
te_c = te_fine[c_ind]

lv = 3*np.sqrt(te_l/te_c)
lw = lv**2

print 'lv = ', lv
print 'lw = ', lw
print 'buffer', l_buffer_x, r_buffer_x
