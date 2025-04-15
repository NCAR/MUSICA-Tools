#! /usr/bin/env python

from netCDF4 import Dataset
from os import getenv

Output_File = getenv('output_file')

F0 = Dataset(Output_File,mode="r+")

print("Grabbing and adding SOAG variables")
SOAG0 = F0.variables['SOAG0'][:,:,:]
SOAG1 = F0.variables['SOAG1'][:,:,:]
SOAG2 = F0.variables['SOAG2'][:,:,:]
SOAG3 = F0.variables['SOAG3'][:,:,:]
SOAG4 = F0.variables['SOAG4'][:,:,:]
SOAG = SOAG0+SOAG1+SOAG2+SOAG3+SOAG4

print("Grabbing and adding soa_a1 variables")
soa1_a1 = F0.variables['soa1_a1'][:,:,:]
soa2_a1 = F0.variables['soa2_a1'][:,:,:]
soa3_a1 = F0.variables['soa3_a1'][:,:,:]
soa4_a1 = F0.variables['soa4_a1'][:,:,:]
soa5_a1 = F0.variables['soa5_a1'][:,:,:]
soa_a1 = soa1_a1+soa2_a1+soa3_a1+soa4_a1+soa5_a1

print("Grabbing and adding soa_a2 variables")
soa1_a2 = F0.variables['soa1_a2'][:,:,:]
soa2_a2 = F0.variables['soa2_a2'][:,:,:]
soa3_a2 = F0.variables['soa3_a2'][:,:,:]
soa4_a2 = F0.variables['soa4_a2'][:,:,:]
soa5_a2 = F0.variables['soa5_a2'][:,:,:]
soa_a2 = soa1_a2+soa2_a2+soa3_a2+soa4_a2+soa5_a2

print("Creating new variables and attributes in file")
vdims=F0.variables['SOAG0'].dimensions
type=F0.variables['SOAG0'].dtype

print("Creating variable SOAG")
F0.createVariable('SOAG',type,vdims)
F0.variables['SOAG'].setncattr('mdims',1)
F0.variables['SOAG'].setncattr('units',"kg/kg")
F0.variables['SOAG'].setncattr('long_name',"SOAG")

vdims=F0.variables['soa1_a1'].dimensions
type=F0.variables['soa1_a1'].dtype

print("Creating variable soa_a1")
F0.createVariable('soa_a1',type,vdims)
F0.variables['soa_a1'].setncattr('mdims',1)
F0.variables['soa_a1'].setncattr('units',"kg/kg")
F0.variables['soa_a1'].setncattr('long_name',"soa_a1")

print("Creating variable soa_a2")
F0.createVariable('soa_a2',type,vdims)
F0.variables['soa_a2'].setncattr('mdims',1)
F0.variables['soa_a2'].setncattr('units',"kg/kg")
F0.variables['soa_a2'].setncattr('long_name',"soa_a2")

print("Writing SOAG to file")
VAR0 = F0.variables['SOAG']
VAR0[:,:,:]=SOAG[:,:,:]

print("Writing soa_a1 to file")
VAR0 = F0.variables['soa_a1']
VAR0[:,:,:]=soa_a1[:,:,:]

print("Writing soa_a2 to file")
VAR0 = F0.variables['soa_a2']
VAR0[:,:,:]=soa_a2[:,:,:]

F0.close()
