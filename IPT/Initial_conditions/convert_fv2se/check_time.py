#! /usr/bin/env python

from netCDF4 import Dataset
from os import getenv

Output_File = getenv('output_file')

F0 = Dataset(Output_File,mode="r+")

print("Grabbing time variables")
time = F0.variables['time'][:]
time_bounds = F0.variables['time_bounds'][:,:]

time=time*0.0
time_bounds=time_bounds*0.0

print("Updateing time and time_bounds")
VAR0 = F0.variables['time']
VAR0[:]=time[:]

VAR1 = F0.variables['time_bounds']
VAR1[:,:]=time_bounds[:,:]

F0.close()
