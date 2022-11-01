#%% User input
# forecasting cycle to be used

# Ian
#cycle = '2022092718'
#storm_num = '09'
#basin = 'al'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2022/sd1059_hurricane_2022.nc'

# Ian
cycle = '2022092718'
storm_num = '09'
basin = 'al'
url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2022/sd1059_hurricane_2022.nc'
#url_saildrone = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/Data/Saildrones/2022/sd1032_hurricane_2022.nc'

# Fiona
#cycle = '2022091806'
#storm_num = '07'
#basin = 'al'
#url_saildrone = scratch_folder + 'Data/Saildrones/2022/sd1031_hurricane_2022.nc'

#exp_names = ['hafsv0p2a_2021rt_natl','HAFSv0p3_pre_baseline','HAFSv0p3_prebs_h3bc']
#exp_labels = ['HAFSv0.2a','HAFSv0.3 prebs control','HAFSv0.3 prebs h3bc']
#exp_colors = ['c','limegreen','indianred']

exp_names = ['HAFSv0p3_test_cpl_bugfix_warm_start','hafsv0p3a_2022rt_natl']
exp_labels = ['HAFSv0.3a_fixed_wind_stresses','HAFSv0.3a_rt']
exp_colors = ['limegreen','c']

#exp_names = ['hafsv0p3a_2022rt_natl']
#exp_labels = ['HAFSv0.3a_rt']
#exp_colors = ['c']

lon_lim = [-75,-55.0]
lat_lim = [20.0,50.0]

home_folder = '/home/Maria.Aristizabal/'
scratch_folder = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'
scratch2_folder = '/scratch2/NCEPDEV/hurricane/noscrub/Maria.Aristizabal/'
abdeck_folder = '/scratch1/NCEPDEV/hwrf/noscrub/input/abdeck/'

folder_exps = [scratch2_folder + exp_names[0] + '/' + cycle + '/' + storm_num + basin[-1] + '/',
               scratch_folder + exp_names[1] + '/' + cycle + '/' + storm_num + basin[-1] + '/']

#folder_exps = [scratch_folder + exp_names[0] + '/' + cycle]

bath_file = scratch_folder +'bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_70.0.nc'

best_track_file = abdeck_folder + 'btk/b' + basin + storm_num + cycle[0:4] + '.dat'
GFS_track_file = abdeck_folder + 'aid/a' + basin + storm_num + cycle[0:4] + '.dat'

# folder utils for Hycom 
#folder_utils4hycom= '/home/Maria.Aristizabal/hafs_graphics/ush/python/ocean/'
folder_utils4hycom= '/home/Maria.Aristizabal/Repos/NCEP_scripts/'
folder_uom= '/home/Maria.Aristizabal/Repos/Upper_ocean_metrics/'
folder_myutils= '/home/Maria.Aristizabal/Utils/'

folder_fig = scratch_folder + 'figures/'

################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import sys
import os
import glob
#import seawater as sw

sys.path.append(folder_utils4hycom)
from utils4HYCOM import readgrids,readdepth,readVar,readBinz

sys.path.append(folder_uom)
from Upper_ocean_metrics import MLD_temp_crit, OHC_from_profile

sys.path.append(folder_myutils)
from my_models_utils import get_storm_track_and_int, get_best_track_and_int,\
                            get_GFS_track_and_int, geo_coord_to_HYCOM_coord,\
                            get_var_from_model_following_trajectory


#plt.switch_backend('agg')

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

################################################################################
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

################################################################################
#%% Time fv3
time_fv3 = [tini + timedelta(hours=int(dt)) for dt in np.arange(0,127,3)]

#################################################################################
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

#################################################################################
#%% Read best track
lon_best_track, lat_best_track, t_best_track, int_best_track, name_storm = get_best_track_and_int(best_track_file)

time_best_track = np.asarray([datetime.strptime(t,'%Y%m%d%H') for t in t_best_track])

#################################################################################
#%% Read GFS track
lon_GFS_track, lat_GFS_track, lead_time_GFS, int_GFS_track, cycle_GFS = get_GFS_track_and_int(GFS_track_file,cycle)
 
#################################################################################
#%% Read Saildrone data
# Sam

url = url_saildrone

gdata = xr.open_dataset(url)#,decode_times=False)

latitude = np.asarray(gdata.latitude)
longitude = np.asarray(gdata.longitude)
time = np.asarray(gdata.time)
dataset_id = gdata.drone_id
temp_air_mean = np.asarray(gdata.TEMP_AIR_MEAN)
rh_mean = np.asarray(gdata.RH_MEAN)
baro_pres_mean = np.asarray(gdata.BARO_PRES_MEAN)
temp_sb37_mean = np.asarray(gdata.TEMP_SBE37_MEAN)
wind_from_mean = np.asarray(gdata.WIND_FROM_MEAN)
wind_speed_mean = np.asarray(gdata.WIND_SPEED_MEAN)
sal_sb37_mean = np.asarray(gdata.SAL_SBE37_MEAN)
#water_current_speed_mean = np.asarray(gdata.WATER_CURRENT_SPEED_MEAN)
#water_current_direccion_mean = np.asarray(gdata.WATER_CURRENT_DIRECTION_MEAN)
wave_dominant_period = np.asarray(gdata.WAVE_DOMINANT_PERIOD)
wave_significant_height = np.asarray(gdata.WAVE_SIGNIFICANT_HEIGHT)

times = np.asarray(gdata.time)
timestamps = mdates.date2num(time)
times = np.asarray(mdates.num2date(timestamps))
oktimeg = np.logical_and(mdates.date2num(time) >= mdates.date2num(tini),\
                         mdates.date2num(time) <= mdates.date2num(tend))

# Fields within time window
latS = latitude[oktimeg]
lonS = longitude[oktimeg]
temp_air = temp_air_mean[oktimeg]
rh = rh_mean[oktimeg]
baro_pres =  baro_pres_mean[oktimeg]
temp_sb37 = temp_sb37_mean[oktimeg]
wind_from = wind_from_mean[oktimeg]
wind_speed = wind_speed_mean[oktimeg]
sal_sb37 = sal_sb37_mean[oktimeg]
#water_speed = water_current_speed_mean[oktimeg]
#water_dir = water_current_direccion_mean[oktimeg]
wave_dom_period = wave_dominant_period[oktimeg]
wave_dom_height = wave_significant_height[oktimeg]
timeS = times[oktimeg]
timestampS = timestamps[oktimeg]

#################################################################################
#%% Loop the experiments

lon_forec_track = np.empty((len(folder_exps),43))
lon_forec_track[:] = np.nan
lat_forec_track = np.empty((len(folder_exps),43))
lat_forec_track[:] = np.nan
lead_time = np.empty((len(folder_exps),43))
lead_time[:] = np.nan
int_track = np.empty((len(folder_exps),43))
int_track[:] = np.nan

#target_timeShyc = np.empty((len(folder_exps),22))
#target_timeShyc[:] = np.nan
target_tempS = np.empty((len(folder_exps),22))
target_tempS[:] = np.nan
target_saltS = np.empty((len(folder_exps),22))
target_saltS[:] = np.nan
#target_timeSfv3 = np.empty((len(folder_exps),43))
#target_timeSfv3 = np.nan
target_press_surf = np.empty((len(folder_exps),43))
target_press_surf[:] = np.nan
target_tmp_surf = np.empty((len(folder_exps),43))
target_tmp_surf[:] = np.nan
target_tmp_2mabove = np.empty((len(folder_exps),43))
target_tmp_2mabove[:] = np.nan
target_rh_2mabove = np.empty((len(folder_exps),43))
target_rh_2mabove[:] = np.nan

for i,folder in enumerate(folder_exps):
    print(folder)
    #%% Get list files
    files_hafs_hycom = sorted(glob.glob(os.path.join(folder,'*3z*.nc')))
    #files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.nc')))
    #if len(files_hafs_fv3) == 0:
    #    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.nc')))
    files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafsprs*.grb2')))
    if len(files_hafs_fv3) == 0:
        files_hafs_fv3 = sorted(glob.glob(os.path.join(folder,'*hafs*grid01*.grb2')))

    #%% Reading HAFS/HYCOM grid
    hafs_hycom_grid = xr.open_dataset(files_hafs_hycom[0],decode_times=False)
    lon_hafs_hycom = np.asarray(hafs_hycom_grid['Longitude'][:])
    lat_hafs_hycom = np.asarray(hafs_hycom_grid['Latitude'][:])
    depth_hafs_hycom = np.asarray(hafs_hycom_grid['Z'][:])

    #%% Reading HAFS/FV3 grid
    #hafs_fv3_grid = xr.open_dataset(files_hafs_fv3[0],decode_times=False)
    #lon_hafs_fv3 = np.asarray(hafs_fv3_grid['longitude'][:])
    #lat_hafs_fv3 = np.asarray(hafs_fv3_grid['latitude'][:])

    fv3 = xr.open_dataset(files_hafs_fv3[0],engine="pynio")
    lon_hafs_fv3 = np.asarray(fv3.lon_0)
    lat_hafs_fv3 = np.asarray(fv3.lat_0)

    #%% Get storm track from trak atcf files
    files = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix')))
    if files:
        file_track = files[0]
    else:
        file_track = sorted(glob.glob(os.path.join(folder,'*trak*.atcfunix.all')))[0]
    okn = get_storm_track_and_int(file_track,storm_num)[0].shape[0]
    lon_forec_track[i,0:okn], lat_forec_track[i,0:okn], lead_time[i,0:okn], int_track[i,0:okn], _ = get_storm_track_and_int(file_track,storm_num)

    #%% Read HAFS/Fv3 time
    '''
    time_fv3 = []
    for n,file in enumerate(files_hafs_fv3):
        print(file)
        FV3 = xr.open_dataset(file)
        t = FV3.variables['time'][:]
        timestamp = mdates.date2num(t)[0]
        time_fv3.append(mdates.num2date(timestamp))

    time_fv3 = np.asarray(time_fv3)
    '''

    #%% Read HAFS/HYCOM time
    time_hycom = []
    timestamp_hycom = []
    for n,file in enumerate(files_hafs_hycom):
        print(file)
        HYCOM = xr.open_dataset(file)
        t = HYCOM['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time_hycom.append(mdates.num2date(timestamp))
        timestamp_hycom.append(timestamp)

    time_hycom = np.asarray(time_hycom)
    timestamp_hycom = np.asarray(timestamp_hycom)

    #################################################################################
    #%% Retrieve HAFS_HYCOM temp. following saildrone trajectory

    # Conversion from glider longitude and latitude to HYCOM convention
    lonS_hyc, latS_hyc = geo_coord_to_HYCOM_coord(lonS,latS)

    files_model = files_hafs_hycom
    time_name = 'MT'
    lat_name = 'Latitude'
    lon_name = 'Longitude'
    depth_level = 0
    timestamp_obss = timestampS
    kwargs = dict(depth_level = 0)
    if np.min(lon_hafs_hycom) < 0:
        lon_obss = lonS
    else: 
        lon_obs = lonS_hyc
    lat_obss = latS_hyc

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    target_timeShyc, target_tempS[i,:] = get_var_from_model_following_trajectory(files_model,'temperature',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

    _, target_saltS[i,:] = get_var_from_model_following_trajectory(files_model,'salinity',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs)

    #################################################################################
    #%% Retrieve HAFS_atm press following saildrone trajectory

    files_model = files_hafs_fv3
    time_name = 'time'
    #lat_name = 'latitude'
    #lon_name = 'longitude'
    depth_level = 0
    if np.min(lon_hafs_fv3) < 0:
        lon_obss = lonS
    else:
        lon_obss = lonS + 360
    lat_obss = latS
    timestamp_obss = timestampS

    oklo = np.isfinite(lon_obss)
    lon_obs = lon_obss[oklo]
    lat_obs = lat_obss[oklo]
    timestamp_obs = timestamp_obss[oklo]

    #target_timeSfv3, target_press = get_var_from_model_following_trajectory(files_model,'PRES_surface',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    #target_timeSfv3, target_tmp_surf = get_var_from_model_following_trajectory(files_model,'TMP_surface',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    #target_timeSfv3, target_tmp_2mabove = get_var_from_model_following_trajectory(files_model,'TMP_2maboveground',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    #target_timeSfv3, target_rh_2mabove = get_var_from_model_following_trajectory(files_model,'RH_2maboveground',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    target_timeSfv3, target_press_surf[i,:] = get_var_from_model_following_trajectory(files_model,'PRES_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    _, target_tmp_surf[i,:] = get_var_from_model_following_trajectory(files_model,'TMP_P0_L1_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    _, target_tmp_2mabove[i,:] = get_var_from_model_following_trajectory(files_model,'TMP_P0_L103_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

    _, target_rh_2mabove[i,:] = get_var_from_model_following_trajectory(files_model,'RH_P0_L103_GLL0',time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs)

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonS, latS,'.-',color='purple',label='Saildrone Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend(loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
# Elsa
if cycle == '2021070100':
    plt.xlim([-85,-40])
    plt.ylim([5,30])
# Ida
if cycle == '2021082800' or cycle =='2021082706':
    plt.xlim([-95,-75])
    plt.ylim([20,40])

#################################################################################
#%% Figure track
lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='silver')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_best_track, lat_best_track,'o-',color='k',label='Best Track')
for i in np.arange(len(exp_names)):
    plt.plot(lon_forec_track[i,::2], lat_forec_track[i,::2],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
plt.plot(lonS, latS,'.-',color='purple',label='Saildrone Track')
#plt.legend(loc='upper right',bbox_to_anchor=[1.3,0.8])
plt.legend() #loc='upper right')
plt.title('Track Forecast ' + storm_num + ' cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.nanmin(lonS)-1,np.nanmax(lonS)+1])
plt.ylim([np.nanmin(latS)-1,np.nanmax(latS)+1])


#################################################################################

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,temp_sb37,'.-',color='blue',label='Saildrone')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeShyc,target_tempS[i,:],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeShyc,target_tempS,'o-',color='c',label='1m below HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
plt.legend()
plt.title('Water Temperature Cycle '+ cycle,fontsize=18)
plt.ylabel('($^oC$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,sal_sb37,'.-',color='blue',label='Saildrone')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeShyc,target_saltS[i,:],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeShyc,target_saltS,'o-',color='c',label='1m below HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
plt.legend()
plt.title('Salinity Cycle '+ cycle,fontsize=18)
plt.ylabel(' ',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,baro_pres,'.-',color='blue',label='Saildrone')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3,target_press_surf[i,:]/100,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeSfv3,target_press/100,'o-',color='c',label='HAFS-HYCOM',markeredgecolor='k',markersize=10,markeredgewidth=2)
plt.legend()
plt.title('PRESS surface Cycle '+ cycle,fontsize=18)
plt.ylabel('(hPa)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,temp_air,'.-',color='blue',label='Saildrone')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3,target_tmp_2mabove[i,:]-273.15,'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeSfv3,target_tmp_2mabove-273.15,'o-',color='c',label='tmp_2maboveground HAFS-HYCOM',markeredgecolor='k',markersize=10)
plt.legend()
plt.title('Air Temperature Cycle '+ cycle,fontsize=18)
plt.ylabel('($^oC$)',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

fig,ax = plt.subplots(figsize=(10, 4))
plt.plot(timeS,rh,'.-',color='b',label='Saildrone')
for i in np.arange(len(exp_names)):
    plt.plot(target_timeSfv3,target_rh_2mabove[i,:],'o-',color=exp_colors[i],markeredgecolor='k',label=exp_labels[i],markersize=7)
#plt.plot(target_timeSfv3,target_rh_2mabove,'o-',color='c',label='2m above HAFS-HYCOM',markeredgecolor='k',markersize=10)
plt.legend()
plt.title('Relative humidity Cycle '+ cycle,fontsize=18)
plt.ylabel('%',fontsize=14)
date_form = DateFormatter("%m-%d")
ax.xaxis.set_major_formatter(date_form)

