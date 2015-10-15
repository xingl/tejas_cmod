import numpy as np


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
