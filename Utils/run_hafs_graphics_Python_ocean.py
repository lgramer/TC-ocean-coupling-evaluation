
import sys
#sys.path.insert(0,'/lfs4/HFIP/hwrfv3/Hyun.Sook.Kim/myconda')
sys.path.append('/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/hafs_graphics/ush/python/ocean/')
sys.path.append('/home/Maria.Aristizabal/Repos/NCEP_scripts/')

from utils4HWRF import readTrack6hrly
from utils import coast180

import xarray as xr
from datetime import datetime, timedelta

import os
import glob

import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

#----------note-----
# The graphdir is defined inside the individual python script.
#graphdir='/mnt/lfs4/HFIP/hwrfv3/Hyun.Sook.Kim/hafsv0p1cpl_ocngraphics'
#-----------------

#######################################################################

model='hafs'
storm='natl'
tcid='09l'
trackon='y'
cycle='2021083000'
COMhafs='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p3_20220521_v0p3a_vida/com/'+cycle +'/'+tcid+'/'

graphdir='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p3_20220521_v0p3a_vida/com/'+cycle +'/'+tcid+'/emc_figures/'

cdir = '/scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'

#######################################################################
'''
model='hafs'
storm='wpac'
tcid='00w'
trackon='y'
cycle='2021062312'

COMhafs='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_wpac/2021062312/00W'

graphdir='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_wpac/2021062312/00W/figures/'

cdir = '/home/Maria.Aristizabal/Repos/hafs_graphics/ush/python/ocean/'
'''
#######################################################################
'''
model='hafs'
storm='natl'
tcid='00l'
trackon='y'
cycle='2019082800'
COMhafs='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2_baseline_Dorian/dorian05l.2019082800_exp/'

graphdir='/home/Maria.Aristizabal/Figures/2019/HAFSv0.2_baseline/Dorian/2019082800/'

cdir = '/home/Maria.Aristizabal/Repos/hafs_graphics/ush/python/ocean/'
'''
#######################################################################
'''
model='hafs'
storm='epac'
tcid='00e'
trackon='y'
cycle='2021062312'
COMhafs='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_epac/2021062312/00E/'

graphdir='/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_epac/2021062312/00E/figures/'

cdir = '/home/Maria.Aristizabal/Repos/hafs_graphics/ush/python/ocean/'
'''
########################################################################
gatcf=glob.glob(COMhafs+'/*.atcfunix.all')

list_storms=[]
list_tcids=[]

#if len(gatcf) >1:
for m,G in enumerate(gatcf):
    tmp=G.partition(COMhafs)
      
    patcf0=tmp[-1].partition('.'+cycle)[0]
    patcf=patcf0.replace('/','')
    list_storms.append(patcf[:-3])
    list_tcids.append(patcf[-3:])
 
#cdir=os.getcwd()
os.chdir(cdir)
#os.chdir(COMhafs)
os.system('module load intel/19.0.5.281')
os.system('module load netcdf/4.7.0')
os.system('module load wgrib2/2.0.8')

#m=0
#cmd='python '+cdir+'/storm_OHC.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir

#os.system(cmd)

#run /scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/hafs_graphics/ush/python/ocean/storm_SST.py hafs natl 09l 2021083000 'y' /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p3_20220521_v0p3a_vida/com/2021083000/09L/ /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p3_20220521_v0p3a_vida/com/2021083000/09L/emc_figures/


#run /scratch1/NCEPDEV/hwrf/save/Maria.Aristizabal/hafs_graphics/ush/python/ocean/SSTnc.py hafs natl 09l 2021083000 'y' /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p3_20220521_v0p3a_vida/com/2021083000/09L/ /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p3_20220521_v0p3a_vida/com/2021083000/09L/emc_figures/

'''
run SSTnc.py hafs wpac 00w 2021062312 'y' /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_wpac/2021062312/00W/ /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_wpac/2021062312/00W/figures/

run SSTnc.py hafs epac 00e 2021062312 'y' /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_epac/2021062312/00E/ /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/hafsv0p2a_2021rt_epac/2021062312/00E/figures/

#run storm_HeatFlux.py hafs natl 00l 2019082800 'y' /scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2_baseline_Dorian/dorian05l.2019082800_exp/ /home/Maria.Aristizabal/Figures/2019/HAFSv0.2_baseline/Dorian/2019082800/
'''

#---- produce basin-scale figures:

'''
cmd='python '+cdir+'SSTnc.py '+model+' '+storm+' '+tcid+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
os.system(cmd) 

cmd='python '+cdir+'/MLDnc.py '+model+' '+storm+' '+tcid+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
os.system(cmd)

cmd='python '+cdir+'/OHCnc.py '+model+' '+storm+' '+tcid+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
os.system(cmd)

cmd='python '+cdir+'/Z20nc.py '+model+' '+storm+' '+tcid+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
os.system(cmd)

#run SSTnc.py model storm tcid cycle trackon /scratch2/AOML/aoml-phod/Maria.Aristizabal/HAFSv0.2_baseline_HYCOM_Dorian/dorian05l.2019082800_exp/ /home/Maria.Aristizabal/Figures/2019/HAFSv0.2_baseline/Dorian/2019082800/

#---- produce storm-scale figures:

for m in range(len(list_storms)):
   cmd='python '+cdir+'/storm_HeatFlux.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_SST.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_MLD.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_OHC.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_Z20.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_tempZ40m.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_tempZ70m.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_tempZ100m.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_WvelZ40m.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_WvelZ70m.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_WvelZ100m.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)

   cmd='python '+cdir+'/storm_HeatFlux.py '+model+' '+list_storms[m]+' '+list_tcids[m]+' '+cycle+' '+trackon+' '+COMhafs+' '+graphdir
   os.system(cmd)
'''
