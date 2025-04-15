#! /usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from os import getenv

def destagger_uv(u_staggered, v_staggered):
    """Destaggers U and V wind components from a staggered grid to cell centers.

    Args:
        u_staggered (numpy.ndarray): U component of wind on the staggered grid.
        v_staggered (numpy.ndarray): V component of wind on the staggered grid.

    Returns:
        tuple: A tuple containing the destaggered U and V components.
    """
    
    print("U size: ",u_staggered.shape)
    print("V size: ",v_staggered.shape)

    #u_destaggered = (u_staggered[:,:,:-1,:] + u_staggered[:,:,1:,:]) / 2
    #u_destaggered = (u_staggered[:,:,:,:-1] + u_staggered[:,:,:,1:]) / 2

    u_destaggered = u_staggered[:,:,:,:]
    u_destaggered = np.append(u_destaggered,u_destaggered[:,:,-2:-1,:]*0.97,axis=2)

    #v_destaggered = (v_staggered[:,:,:-1,:] + v_staggered[:,:,1:,:]) / 2
    #v_destaggered = (v_staggered[:,:,:,:-1] + v_staggered[:,:,:,1:]) / 2    

    v_destaggered = v_staggered

    print("U size new: ",u_destaggered.shape)
    print("V size new: ",v_destaggered.shape)

    return u_destaggered, v_destaggered


#Input_File = getenv('input_file')
Output_File = getenv('output_file')

with Dataset(Output_File, "r+") as ncin:
    US = ncin.variables['US'][:,:,:,:]
    VS = ncin.variables['VS'][:,:,:,:]
    U,V = destagger_uv(US,VS)
        
    vars=['US','VS']
    vars_new=['U','V']
    for var,var_new in zip(vars,vars_new):
        print("Check for: ",var,var_new)
        type=ncin.variables[var].dtype #note that case makes a difference
        vdims=ncin.variables['T'].dimensions
        ncin.createVariable(var_new,type,vdims)
        varatts = ncin.variables[var].ncattrs()
        for att in varatts:
            if att != 'long_name' :
                val = ncin.variables[var].getncattr(att)
                ncin.variables[var_new].setncattr(att,val)
            else:
                if var == 'US':
                    ncin.variables[var_new].setncattr(att,'Zonal wind')
                if var == 'VS':
                    ncin.variables[var_new].setncattr(att,'Meridional wind')

    U0 = ncin.variables['U']
    U0[:,:,:,:]=U[:,:,:,:]
        
    V0 = ncin.variables['V']
    V0[:,:,:,:]=V[:,:,:,:]


