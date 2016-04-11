#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv,exit,stdout
import optparse as op
from ParIO import *
import os
import sys
import optparse as op
from subprocess import call
from read_write_geometry import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

parser = op.OptionParser()
options,args = parser.parse_args()
geom_file = args[0]
kymin = float(args[1])

gpars,geometry = read_geometry_global(geom_file)

nz,nx=np.shape(geometry['gxx'])
#geom_list = ['gyy']
#for i in geom_list :
#    geom = geometry[i]
    #z_ind = np.argmin(abs(geom[:,1]))
    #print z_ind
#    geom_obmp=geom[32,:]
#    plt.plot(geom_obmp)
#    plt.title(i)
#    plt.show()

q =  gpars['q0']
Lref = gpars['Lref']
cy = geometry ['C_y'][1]
gyy = geometry ['gyy'][32][1]
gyz = geometry ['gyz'][32][1]
gzz = geometry ['gzz'][32][1]
print 'q, Cy, gyy, gyz, gzz'
print q, cy, gyy, gyz, gzz

#kymin =   0.5000E-01
rhostar  =   0.16839125E-02
print Lref*rhostar
k_theta = q*cy/np.sqrt(q**2*cy**2*gyy+2*q*cy*gyz+gzz)*kymin/Lref/rhostar
print k_theta
