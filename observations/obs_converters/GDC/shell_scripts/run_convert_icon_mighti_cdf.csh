#!/bin/csh

#@ i = 22
#@ j = 31
#while ($i <= $j) 
#   set obstime = `printf "2009-03-%02d_23:00:00" $i`
#   echo $obstime
#   set datestr  = " datestr = "$obstime
#   cp input.nml input.nml.original
#   sed -i  's/'"`grep 'datestr' input.nml.original`"'/'"$datestr"'/g' input.nml # >! input.nml
# 
#   ./convert_cosmic_gps_ncf
#
#   @ i += 1
#end
#rm convert.nml

set yyyy = 2020
set mm = 04

@ iday = 1
@ jday = 30

while ($iday <= $jday)
  @ i = 00
  @ j = 23
  while ($i <= $j) 

    #input files:
    set date = `printf "%04d-%02d-%02d" $yyyy $mm $iday`
    ls /glade/scratch/nickp/tmp/icon/ICON_L2-2_MIGHTI_Vector-Wind-Green_${date}_v04*NC >! filelist.txt
    ls /glade/scratch/nickp/tmp/icon/ICON_L2-2_MIGHTI_Vector-Wind-Red_${date}_v04*NC >> filelist.txt

    set obstime = `printf "%04d-%02d-%02d_%02d:00:00" $yyyy $mm $iday $i`
    set outputfile = `printf "obs_seq.icon.%04d%02d%02d%02d" $yyyy $mm $iday $i`
    echo $obstime
    set datestr          = " datestr         = "$obstime
    set mighti_out_file  = " mighti_out_file = "$outputfile
    echo $mighti_out_file
    sed -e  's/'"`grep 'datestr'         input.nml.original`"'/'"$datestr"'/' \
        -e  's/'"`grep 'mighti_out_file' input.nml.original`"'/'"$mighti_out_file"'/' \
    input.nml.original >!  input.nml
    ./convert_icon_mighti_ncf 
  
    mv $outputfile /glade/scratch/nickp/tmp/icon/obs_seq/
 
    @ i += 1
  end
  @ iday += 1
end
#mv obs_seq.out obs_seq.out.23






