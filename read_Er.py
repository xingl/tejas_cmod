import numpy as np

def read_Er(Er_file, shift_Er, Er_shift):
 
    Er_data = np.genfromtxt(Er_file)
    if shift_Er:
       Er_data[:,0] = Er_data[:,0]+Er_shift
    
    psip_norm_Er = Er_data[:,0]
    Er = Er_data[:,1]
    Er_error = Er_data[:,2]

    return psip_norm_Er, Er, Er_error
       
