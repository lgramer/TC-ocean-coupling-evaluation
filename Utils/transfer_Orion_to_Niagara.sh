#!/bin/bash

module load globus-cli
module load hpss
/home/Maria.Aristizabal/load_miniconda.sh

export ORION=84bad22e-cb80-11ea-9a44-0255d23c44ef
export NIAGARA=21467dd0-afd6-11ea-8f12-0a21f750d19b

globus login
globus endpoint activate --web --no-browser $ORION
globus endpoint activate --web --no-browser $NIAGARA

globus ls $ORION\:/work/noaa/hwrf/noscrub/maristiz/hafsarch/hafsv0p2a_phase3/
globus ls $NIAGARA\:/collab1/data/Maria.Aristizabal/hafsv0p2a_phase3/

python transfer_Orion_to_Niagara.py
