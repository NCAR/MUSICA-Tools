#! /usr/bin/env python

from netCDF4 import Dataset
from os import getenv
import numpy as np

Output_File = getenv('output_file')

F0 = Dataset(Output_File,mode="r+")

with Dataset(Output_File, 'r') as ds:
    for var_name in ds.variables:
        var_data = ds.variables[var_name][:]

        try:
            if np.isnan(var_data).any():
                print(f"NaNs found in variable: {var_name}")
            if np.isinf(var_data).any():
                print(f"Infs found in variable: {var_name}")
        except: 
            print(var_name," data type not supported")

F0.close()

