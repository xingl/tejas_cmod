#!/usr/bin/env python

import numpy as np
import re
import scipy as sc
import matplotlib.pyplot as plt
from scipy import interpolate

def iterdb_header(quant_str,units,quant_len,shot_num):
    header=' '+shot_num+' 2 0 6             ;-SHOT #- F(X) DATA -gyro2iterdb-\n'
    header=header+'                              ;-SHOT DATE-  UFILES ASCII FILE SYSTEM\n'
    header=header+'   0                          ;-NUMBER OF ASSOCIATED SCALAR QUANTITIES\n'
    header=header+' RHOTOR              -        ;-INDEPENDENT VARIABLE LABEL: X-\n'
    header=header+' TIME                SECONDS  ;-INDEPENDENT VARIABLE LABEL: Y-\n'
    header=header+' '+quant_str+'              '+units+'        ;-DEPENDENT VARIABLE LABEL\n'
    header=header+' 3                            ;-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM\n'
    header=header+'         '+str(quant_len)+'                   ;-# OF X PTS- \n'
    header=header+'          1                   ;-# OF Y PTS-  X,Y,F(X,Y) DATA FOLLOW:\n'
    return header 

def iterdb_write_quant(fileid,quant_arr):
    fileid.write(' ')
    for i in range(len(quant_arr)):
        out_str='%-12e' % quant_arr[i]
        if quant_arr[i] > 0.0:
            out_str=' '+out_str
        fileid.write(out_str)
        if (i+1)%6==0 and i != len(quant_arr)-1:
            fileid.write('\n')
            fileid.write(' ')
    fileid.write('\n')

def output_iterdb(rhot,ne,te,ni,ti,file_base,shot_num,time_string,vrot=np.array(0),nimp=np.array(0)):
    # takes te, ti in the unit of ev
    # takes ne, ni in the unit of m^(-3)
    transition=';----END-OF-DATA-----------------COMMENTS:-----------\n'
    transition=transition+'********************************************************************************\n'
    transition=transition+'********************************************************************************\n'
    idbf=open(file_base+'.iterdb','w')
    header=';Created with script profiles_pol2tor.py for GENE input\n' 
    idbf.write(header)
    header=';----END-OF-ORIGINAL-HEADER------COMMENTS:-----------\n'
    idbf.write(header)

    header=iterdb_header('TE','ev',len(rhot),shot_num)
    idbf.write(header)
    iterdb_write_quant(idbf,rhot)
    idbf.write('  '+time_string+'\n')      
    iterdb_write_quant(idbf,te)
    idbf.write(transition)

    header=iterdb_header('TI','ev',len(rhot),shot_num)
    idbf.write(header)
    iterdb_write_quant(idbf,rhot)
    idbf.write('  '+time_string+'\n')      
    iterdb_write_quant(idbf,ti)
    idbf.write(transition)

    header=iterdb_header('NE','m-3',len(rhot),shot_num)
    idbf.write(header)
    iterdb_write_quant(idbf,rhot)
    idbf.write('  '+time_string+'\n')      
    iterdb_write_quant(idbf,ne)
    idbf.write(transition)

    header=iterdb_header('NM1','m-3',len(rhot),shot_num)
    idbf.write(header)
    iterdb_write_quant(idbf,rhot)
    idbf.write('  '+time_string+'\n')      
    iterdb_write_quant(idbf,ni)
    idbf.write(transition)

    #Include impurity density
    if nimp.any():
        header=iterdb_header('NM2','m-3',len(rhot),shot_num)
        idbf.write(header)
        iterdb_write_quant(idbf,rhot)
        idbf.write('  '+time_string+'\n')      
        iterdb_write_quant(idbf,nimp)
        idbf.write(transition)

    if vrot.any():
        header=iterdb_header('VROT','rad/s',len(rhot),shot_num)
        idbf.write(header)
        iterdb_write_quant(idbf,rhot)
        idbf.write('  '+time_string+'\n')      
        iterdb_write_quant(idbf,vrot)
        idbf.write(transition)

    idbf.close()


