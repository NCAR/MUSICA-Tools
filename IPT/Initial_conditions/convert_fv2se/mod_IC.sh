#!/bin/bash -l

#####################################################################
#
#Script Name : mod_IC.sh
#
#Description : This script converts an FV IC file to a SE IC file.
#              conda, NPL, and NCO are required. 
#
# Input Args : n/a
#
# Author     : Shawn Honomichl
#
# Modification History :
#              11 December, 2024 - Primary script written by
#              Shawn Honomichl.
#
##########################################################################
#                          EDITABLE REGION                               #
##########################################################################

#note that the input file for this script should be an Finite Volume (FV) IC file
input_file=/glade/campaign/acom/acom-climate/UTLS/shawnh/1980-01-01-00000/f.e21.FWHISTBgcCrop.f09_f09_mg17.CMIP6-AMIP-WACCM.001.cam.i.1980-01-01-00000.nc

#set a full path output filename - temporary template file will also get written to this directory
output_file=/glade/campaign/acom/acom-climate/UTLS/shawnh/1980-01-01-00000/WACCMFV2CAMSE.V2.cam.i.1980-01-01-00000.nc

#set a spectral element (SE) IC file to use to regrid the FV IC file to an SE grid.
#template_file=/glade/campaign/acom/acom-climate/UTLS/shawnh/archive/f.e30_beta02.FCts4MTHIST.ne30.104/atm/hist/f.e30_beta02.FCts4MTHIST.ne30.104.cam.i.1981-01-01-00000.nc
template_file=/glade/u/home/shawnh/Mod_IC/GRID_107L_WACCM.nc

#set the remapping grid to convert the FV grid to SE
map_grid_file=/glade/work/aherring/grids/tempestwgts/f09_2_ne30_fv2se_stt.nc

##########################################################################
#                      END OF EDITABLE REGION                            #
##########################################################################

#load modules 
module load conda/latest #load up the latest python conda environment
conda activate npl #load the NCAR package library
module load nco #load nco for regridding

#export input_file and output_file for python use
export input_file=$input_file
export output_file=$output_file

#horizontal grid interpolation from FV -> SE
echo +++++++++++++++++++++++++++++++++++++++++++
echo converting input file from FV to ne30 grid
echo "in new output file "$output_file
echo +++++++++++++++++++++++++++++++++++++++++++
ncremap --add_fll -m $map_grid_file $input_file $output_file
echo done converting input file
echo " "

#rename variables to match the format of the ne30 format
echo ++++++++++++++++++++++++++++++++++
echo Renaming ncol, lat, lon, and area
echo ++++++++++++++++++++++++++++++++++
ncrename -v lat,lat_d -v lon,lon_d -v area,area_d -d ncol,ncol_d -O $output_file
echo Done renaming variables
echo " "

#re-introduce the lat_d and lon_d fields from the template file
#since ncremap doesn't exacly map the lat_d field correctly.
#Otherwise CESM will spit out an error:
#ncdata file latitudes not in correct column order
echo ++++++++++++++++++++++++++++++++++++++++++++++
echo Re-copying lat_d and lon_d from template file
echo to prevent CESM lat_d column order error
echo ++++++++++++++++++++++++++++++++++++++++++++++
ncks --chk_nan -A -v lat_d,lon_d $template_file $output_file
echo Done re-copying lat_d and lon_d from template file
echo " "

#copy over the template file, and add P0 to it.
echo +++++++++++++++++++++++++++++++++++++++++
echo copying over template file and adding P0
echo for vertical interpolation step
echo +++++++++++++++++++++++++++++++++++++++++
cp -f $template_file $output_file.template.nc #copy file over and rename it
template_file=$output_file.template.nc
ncks --chk_nan -A -v P0 $input_file $template_file #add P0 to template file
echo done copying template file and adding P0
echo " "

#vertical interpolation to match
#the template file's vertical levels
echo ++++++++++++++++++++++++++++++++
echo converting vertical grid 
echo ++++++++++++++++++++++++++++++++
ncremap --vrt=$template_file $output_file $output_file.tmp.nc
rm $output_file
mv $output_file.tmp.nc $output_file
echo done converting vertical grid of input file
echo " "

#temporarily convert the file to NetCDF4-Classic format so that the remaining
#operations can be done faster and more efficiently.  Otherwise
#the python program that adds variables will go painfully slow.
echo ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
echo Temp converting IC file to NetCDF4 Classic so that python runs faster
echo ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ncks --chk_nan --fl_fmt=netcdf4_classic -O --no_alphabetize $output_file $output_file
echo done converting to NetCDF4
echo " "

#Calculate and add in the SOAG, soa_a1, and soa_a2 variables
#to the output file.
echo +++++++++++++++++++++++++++++++++++++++++++++++
echo adding SOAG, soa_a1, and soa_a2 to output file
echo +++++++++++++++++++++++++++++++++++++++++++++++
python3 mod_vars.py
echo done adding SOAG, soa_a1, and soa_a2 to output file
echo " "

#Add in U and V winds 
echo ++++++++++++++++++++++++++++++++++++
echo Adding U/V winds from template file
echo ++++++++++++++++++++++++++++++++++++
ncks --chk_nan -A -v U,V $template_file $output_file
rm $template_file
echo done adding U/V winds
echo " "

#Do a time check to see if the time and time bounds variables are off.
#if so then correct the time
echo +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
echo Checking consistency of date/time variables in output file
echo +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
python3 check_time.py
echo date/time variable check done
echo " "

#perform a NaN/Inf check on the new output file to make sure the file
#is error-free
echo +++++++++++++++++++++++++++++++++++++++
echo Checking output file for NaNs and Infs
echo +++++++++++++++++++++++++++++++++++++++
python3 check_for_naninf.py
echo NaN/Inf check complete
echo " "

#convert back to the original file format when done
#with python operations
echo +++++++++++++++++++++++++++++++++++++++++
echo Converting IC file back to 64 bit offset
echo +++++++++++++++++++++++++++++++++++++++++
ncks --chk_nan --fl_fmt=64bit_offset -O --no_alphabetize $output_file $output_file
echo done converting to 64 bit offset
echo " "
