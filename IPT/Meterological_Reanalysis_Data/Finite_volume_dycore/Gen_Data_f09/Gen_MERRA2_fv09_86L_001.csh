#!/bin/csh
#
#!/bin/tcsh
#PBS -N Gen_MERRA2_86L
#PBS -A P93300043
#PBS -l select=2:ncpus=4:mpiprocs=4:mem=20GB
#PBS -l walltime=12:00:00
#PBS -q casper
#PBS -o logs.gen_MERRA2_86L.out
#PBS -e err.gen_MERRA2_86L.out
#
#---------------------------------------
# Script to Generate MERRA2 Nudging data
#---------------------------------------
module load ncl

#===================================================================
# CONFIGRATION SECTION:
#===================================================================

# Set a REFERENCE (Starting) Date and the number of days to process
#-------------------------------------------------------------------
set RUNNUM   = 01 
set REF_DATE = '20190101'
set NUM_DAYS = 1

# Set INPUT/OUTPUT/TMP directories, 
#-----------------------------------
set NAMELIST = './Config/Config_makeIC-'$RUNNUM'.nl'
set MYLOGDIR  = '/glade/derecho/scratch/shawnh/GEOS5_frcst_data/LOG/LOG_001.'$RUNNUM'/'
set MYTMPDIR  = '/glade/derecho/scratch/shawnh/GEOS5_frcst_data/TMP/TMP_001.'$RUNNUM'/'
set MYOUTDIR  = '/glade/derecho/scratch/sahwnh/'

set INPUTDIR  = '/glade/campaign/collections/rda/data/d313003/orig_res/'
#incomplete:
#set INPUTDIR  = '/glade/campaign/collections/rda/data/ds313.3/0.9x1.25/'


# Set Clean option
#------------------------
set TMP_clean   = 'True'

# Set Processing options
#------------------------
set CASE                   = 'MERRA2_f09_L86'
set DYCORE                 = 'fv'
set PRECISION              = 'float'
set VORT_DIV_TO_UV         = 'False'
set SST_MASK               = 'False'
set ICE_MASK               = 'False'
set OUTPUT_PHIS            = 'True'
set REGRID_ALL             = 'False'
set ADJUST_STATE_FROM_TOPO = 'True'

# Set files containig OUTPUT Grid structure and topography
#---------------------------------------------------------
set fname_grid_info    = '/glade/work/emmons/ic/template_f09_86L.nc'

set fname_phis_output  = '/glade/campaign/cesm/cesmdata/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc'
set ftype_phis_output  = 'FV_TOPOGRAPHY'

# Set INPUT filename format, type of file, and number of time records per file
#------------------------------------------------------------------------------
set fname   = ( none none none none none none )
set ftype   = ( none none none none none none )
set ftime   = ( none none none none none none )

# MERRA (every 3 hours)
set hoursec = ( 00000 10800 21600 32400 43200 54000 64800 75600)
set hourstr = (   00    03    06    09    12    15    18     21   )
set fname[1] = "YYYY/MERRA2_orig_res_YYYYMMDD.nc"
#set fname[1] = "YYYY/MERRA2_0.9x1.25_YYYYMMDD.nc"
set ftype[1] = "MERRA2"
set ftime[1] = "8X"

#GEOS
#set hoursec = ( 00000 21600 43200 64800)
#set hourstr = (   00    06    12    18  )
#set fname[1] = "YYYYMMDD/GEOS5_orig_res_YYYYMMDD+YYYYMMDD.nc"
#set fname[1] = "YYYY/GEOS5_orig_res_YYYYMMDD.nc"

# Set the OUTPUT fileds and the index of the file containing them
#------------------------------------------------------------------
set fields       = ( U V US VS T Q  PS )
set field_findex = ( 1 1 1 1 1 1  1  )

# Set the index of the file containing the INPUT topography
#------------------------------------------------------------------
set phis_findex  = 1


#===================================================================
# PROCESSING SECTION:
#   Loop over range of times, 
#   generate inpute namelist for each time, 
#   then generate the desired Nudging data
#===================================================================

# Construct the fields string
#-----------------------------
set fieldstr  = $fields[1]
set jj        = $field_findex[1]
@ jj = $jj - 1 
set findexstr = $jj
set ii = 1
set fnum = $#fields
while( $ii <  $fnum )
  @ ii = $ii + 1
  set jj = $field_findex[$ii]
  @ jj = $jj - 1
  set fieldstr  = `echo $fieldstr','$fields[$ii]`
  set findexstr = `echo $findexstr','$jj`
end

set jj = $phis_findex
@ jj = $jj - 1 
set phisstr = $jj

# Build NCL fortran library
#---------------------------
if( ! -e MAKEIC.so ) then
#/usr/local/bin/WRAPIT         MAKEIC.stub MAKEIC.f90
#/contrib/ncl-6.1.0/bin/WRAPIT MAKEIC.stub MAKEIC.f90
/glade/u/apps/opt/ncl/6.6.2/bin/WRAPIT MAKEIC.stub MAKEIC.f90
endif

#
# Loop over range of times to process
#  For the first time true, force 
#  generation of ESMF weight datasets
#=====================================
set NDAYS = 0
while ( $NDAYS < $NUM_DAYS )

  # Set current time values
  #------------------------
  set Yearstr = "`date --date=$REF_DATE+${NDAYS}day +%Y`"
  set Monstr  = "`date --date=$REF_DATE+${NDAYS}day +%m`"
  set Daystr  = "`date --date=$REF_DATE+${NDAYS}day +%d`"
  echo $Yearstr ' / ' $Monstr ' / ' $Daystr

  # Loop over the hourly values (4X daily)
  #========================================
  foreach hnum ( 1 2 3 4 5 6 7 8 )
    
    # Set Values dependend upon $hnum, clean 
    # up TMP files at the end of each DAY
    #----------------------------------------------------------------------
    set datestr = $Yearstr$Monstr$Daystr$hoursec[$hnum]
    set LOGFILE = $MYLOGDIR'/LogNCL.'$Yearstr$Monstr$Daystr$hoursec[$hnum]
    if( $hnum == 8 ) then
      set TMP_clean   = 'True'
    else
      set TMP_clean   = 'False'
    endif

    # Create Namelist file for current time
    #----------------------------------------
    echo '                                       '             > $NAMELIST
    echo '! Generated Namelist for makeIC.ncl '            >> $NAMELIST
    echo '!-------------------                   '            >> $NAMELIST
    echo '&makeic_nl'                                         >> $NAMELIST
    echo ' MYTMPDIR="'$MYTMPDIR'"'                            >> $NAMELIST
    echo ' MYOUTDIR="'$MYOUTDIR'"'                            >> $NAMELIST
    echo ' INPUTDIR="'$INPUTDIR'"'                            >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' TMP_clean  ="'$TMP_clean'"'                        >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' REF_DATE              ="'$REF_DATE'"'              >> $NAMELIST
    echo ' CASE                  ="'$CASE'"'                  >> $NAMELIST
    echo ' DYCORE                ="'$DYCORE'"'                >> $NAMELIST
    echo ' PRECISION             ="'$PRECISION'"'             >> $NAMELIST
    echo ' VORT_DIV_TO_UV        ="'$VORT_DIV_TO_UV'"'        >> $NAMELIST
    echo ' SST_MASK              ="'$SST_MASK'"'              >> $NAMELIST
    echo ' ICE_MASK              ="'$ICE_MASK'"'              >> $NAMELIST
    echo ' OUTPUT_PHIS           ="'$OUTPUT_PHIS'"'           >> $NAMELIST
    echo ' REGRID_ALL            ="'$REGRID_ALL'"'            >> $NAMELIST
    echo ' ADJUST_STATE_FROM_TOPO="'$ADJUST_STATE_FROM_TOPO'"'>> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' fname_phis_output     ="'$fname_phis_output'"'     >> $NAMELIST
    echo ' ftype_phis_output     ="'$ftype_phis_output'"'     >> $NAMELIST
    echo ' fname_grid_info       ="'$fname_grid_info'"'       >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' fields       ="'$fieldstr'"'                       >> $NAMELIST
    echo ' source_files ="'$findexstr'"'                      >> $NAMELIST
    echo ' fname_phis_in='$phisstr                            >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      set file = $fname[$ii]
      set file = `echo $file | sed 's/YYYY/'$Yearstr'/g'`
      set file = `echo $file | sed 's/MM/'$Monstr'/g'`
      set file = `echo $file | sed 's/DD/'$Daystr'/g'`
      if( $ftime[$ii] == "4X" ) then
        set file = `echo $file | sed 's/HH/'$hourstr[1]'/'`
        set file = `echo $file | sed 's/HH/'$hourstr[4]'/'`
      else
        set file = `echo $file | sed 's/HH/'$hourstr[$hnum]'/g'`
      endif
      echo " fname$jj="'"'$file'"'                      >> $NAMELIST
    end
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      echo " ftype$jj="'"'$ftype[$ii]'"'                      >> $NAMELIST
    end
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      if( $fname[$ii] == "none" ) then
        echo " fdate$jj="'"-1"'                               >> $NAMELIST
      else
        echo " fdate$jj="'"'$datestr'"'                       >> $NAMELIST
      endif
    end
    echo "/"                                                  >> $NAMELIST

    # Execute NCL program to process the data
    #----------------------------------------
    cat $NAMELIST                     >& $LOGFILE.cfg
    ncl makeIC.ncl NameNumber=$RUNNUM >& $LOGFILE

  # End Loop over the hourly values (4X daily)
  #===========================================
  end # foreach hnum ( 1 2 3 4 )

# End loop over time
#=======================
@ NDAYS = $NDAYS + 1
end # while ( $NDAYS < $NUM_DAYS )

