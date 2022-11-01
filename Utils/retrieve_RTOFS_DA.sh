#!/bin/sh

module load hpss

#date -d "Sun Sep 11 07:59:16 IST 2012+10 days"
date_ini=$(date -d "20201002" +"%Y%m%d")
echo $date_ini

days='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25' #26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60'
for i in $days
do
    cycle=$(date -d "$date_ini +$i days" +"%Y%m%d")
    echo $i
    echo $cycle

    RTOFS_DA_dir='rtofs.'$cycle
    echo $RTOFS_DA_dir
    
    cd /scratch2/AOML/aoml-phod/Maria.Aristizabal/RTOFS-DA/
    mkdir ${RTOFS_DA_dir}
    cd ${RTOFS_DA_dir}

    hsi get /NCEPDEV/emc-ocean/5year/emc.ncodapa/parallel/${RTOFS_DA_dir}/rtofs.ncgrb.tar 

    tar -xf rtofs.ncgrb.tar -v --wildcards 'rtofs_glo_3dz_f006_6hrly_hvr_US_east.nc'
    tar -xf rtofs.ncgrb.tar -v --wildcards 'rtofs_glo_3dz_f012_6hrly_hvr_US_east.nc'
    tar -xf rtofs.ncgrb.tar -v --wildcards 'rtofs_glo_3dz_f018_6hrly_hvr_US_east.nc'
    tar -xf rtofs.ncgrb.tar -v --wildcards 'rtofs_glo_3dz_f024_6hrly_hvr_US_east.nc'
 
    rm rtofs.ncgrb.tar 

done

