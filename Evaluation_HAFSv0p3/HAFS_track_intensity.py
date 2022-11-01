#%% User input

# forecasting cycle to be used
#cycle = '2021070100'
#storm_num = '05'

cycle = '2022092718'
storm_num = '09'
basin = 'al'

exp_names = ['HAFSv0p3_test_cpl_bugfix_warm_start','hafsv0p3a_2022rt_natl']
exp_labels = ['HAFSv0.3a_fixed_wind_stresses','HAFSv0.3a_rt']
exp_colors = ['limegreen','c']

#exp_names = ['hafsv0p2a_2021rt_natl','HAFSv0p3_baseline','HAFSv0p3_HFAB_HF3A']
#exp_labels = ['HAFSv0.2a','H3BL','HF3A']
#exp_colors = ['limegreen','indianred','orange']

lon_lim = [-98.5,-70.0]
lat_lim = [15.0,40.0]

home_folder = '/home/Maria.Aristizabal/'
scratch1_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch2_folder = '/scratch2/NCEPDEV/hurricane/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch2_folder + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/',
               scratch1_folder + exp_names[1] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

bath_file = scratch1_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'


bath_file = scratch1_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = scratch1_folder + 'adeck/a' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = scratch1_folder + 'figures/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            Haversine, get_glider_transect_from_HAFS_HYCOM,\
                            figure_transect_temp, glider_data_vector_to_array,\
                            grid_glider_data


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

###################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

###################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

###################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

###################################################################
#%% Loop the experiments
lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan
for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

###################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
#%% Read GFS track
#lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)
 
##########################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)
okt = np.logical_and(time_best_track >= time_fv3[0],time_best_track <= time_fv3[-1])

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
for i in np.arange(len(exp_names)): 
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.legend(loc='upper right')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
#plt.legend()
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.min(lon_best_track[okt])-5,np.max(lon_best_track[okt])+5])
plt.ylim([np.min(lat_best_track[okt])-5,np.max(lat_best_track[okt])+5])

########################################################################
#%% Figure intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)
tt_fv3 = np.asarray([t.replace(tzinfo=None) for t in time_fv3])
oktt = np.logical_and(tt_fv3[::2] >= time_best_track[0],tt_fv3[::2] <= time_best_track[-1])

fig,ax1 = plt.subplots(figsize=(7, 4))
#plt.plot(lead_time[i,::2][oktt],int_best_track[okt],'o-k',label='Best')
#plt.plot(lead_time[i,::2],int_best_track[okt],'o-k',label='Best')
for ii,i in enumerate(np.arange(len(exp_names))): 
    plt.plot(lead_time[i,::2],int_track[i,::2],'o-',color=exp_colors[i],label=exp_labels[i],markeredgecolor='k',markersize=7)

plt.plot(lead_time[i,::2][oktt],int_best_track[okt],'o-k',label='Best')

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
#ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.legend(loc='upper right',fontsize=14)
#plt.legend(loc='lower right',fontsize=14)
plt.ylim([10,140])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast Cycle '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([10,140])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)

#################################################################################
