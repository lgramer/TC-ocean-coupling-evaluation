module load globus-cli
globus login

export ORION=84bad22e-cb80-11ea-9a44-0255d23c44ef
export NIAGARA=21467dd0-afd6-11ea-8f12-0a21f750d19b
export HERA=82109590-c090-11ea-bef9-0e716405a293
export JET=34ea8506-1882-11eb-81b5-0e2f230cc907

globus endpoint show $ORION | grep Activated
globus endpoint show $NIAGARA | grep Activated
globus endpoint show $HERA | grep Activated
globus endpoint show $JET | grep Activated

globus endpoint activate --web --no-browser $ORION
globus endpoint activate --web --no-browser $NIAGARA
globus endpoint activate --web --no-browser $HERA
globus endpoint activate --web --no-browser $JET

#globus ls $ORION\:/work/noaa/hwrf/noscrub/input/abdeck/aid/
#globus ls $NIAGARA\:/collab1/data/Maria.Aristizabal/
#globus ls $JET\:/mnt/lfs4/HFIP/hwrf-data/hwrf-input/abdeck/btk/

globus transfer $ORION\:/work/noaa/hwrf/noscrub/input/abdeck/aid/aal132020.dat $HERA\:/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/adeck/aal132020.dat
#globus transfer --notify off $JET\:/mnt/lfs4/HFIP/hwrf-data/hwrf-input/abdeck/btk/bwp992020.dat $HERA\:/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Best_track_data/bwp992020.dat


#export prefix=bal
#export year=2020
#for n in  00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
#do 
#globus transfer --notify off $JET\:/mnt/lfs4/HFIP/hwrf-data/hwrf-input/abdeck/btk/${prefix}${n}${year}.dat $HERA\:/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Best_track_data/${prefix}${n}${year}.dat
#done

