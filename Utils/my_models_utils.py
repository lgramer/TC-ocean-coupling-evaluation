
import numpy as np
import matplotlib.dates as mdates
import glob
import os
import xarray as xr
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

################################################################################
#%% Get storm track from trak atcf files
def get_storm_track_and_int(file_track,storm_num):
 
    ff = open(file_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    lead_time = []
    intt = []
    rmww = []
    for l in f:
        if l.split(',')[1].strip() == storm_num:
            lat = float(l.split(',')[6][0:4])/10
            if l.split(',')[6][4] == 'N':
                lat = lat
            else:
                lat = -lat
            lon = float(l.split(',')[7][0:5])/10
            if l.split(',')[7][4] == 'E':
                lon = lon
            else:
                lon = -lon
            latt.append(lat)
            lont.append(lon)
            lead_time.append(int(l.split(',')[5][1:4]))
            intt.append(float(l.split(',')[8]))
            rmww.append(float(l.split(',')[19]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    intt = np.asarray(intt)
    rmww = np.asarray(rmww)
    lead_time, ind = np.unique(lead_time,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]
    int_track = intt[ind]
    rmw_track = rmww[ind]

    return lon_track, lat_track, lead_time, int_track, rmw_track

#################################################################################%% Get high resolution storm track from track_d03.patcf files
def get_storm_track_and_int_high_resolution(file_track):

    ff = open(file_track,'r')
    f = ff.readlines()

    lat = []
    lon = []
    lead_time = []
    w10 = []
    pmin = []
    RMW = []
    for l in f:
        latt = float(l.split(',')[3].split('=')[-1][:-1])
        if l.split(',')[3].split('=')[-1][-1] == 'N':
            latt = latt
        else:
            latt = -latt
        lonn = float(l.split(',')[4].split('=')[-1][:-1])
        if l.split(',')[4].split('=')[-1][-1] == 'E':
            lonn = lonn
        else:
            lonn = -lonn
        lat.append(latt)
        lon.append(lonn)
        lead_time.append(float(l.split(',')[0])) # seconds from beg. of cycle
        w10.append(float(l.split(',')[1].split('=')[-1][:-2])*0.514) # m/s
        pmin.append(float(l.split(',')[2].split('=')[-1][:-4])) # mbar
        RMW.append(float(l.split(',')[5].split('=')[-1][:-5])*1.61) # km

    lat = np.asarray(lat)
    lon = np.asarray(lon)
    lead_time = np.asarray(lead_time)
    w10 = np.asarray(w10)
    pmin = np.asarray(pmin)
    RMW = np.asarray(RMW)
    #lead_time, ind = np.unique(lead_time,return_index=True)
    #lat_track = latt[ind]
    #lon_track = lont[ind]
    #int_track = intt[ind]

    return lon, lat, lead_time, w10, pmin, RMW

################################################################################
def get_best_track_and_int(file_best_track):

    ff = open(file_best_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    time = []
    intt = []
    for l in f:
        lat = float(l.split(',')[6][0:4])/10
        if l.split(',')[6][-1] == 'N':
            lat = lat
        else:
            lat = -lat
        lon = float(l.split(',')[7][0:5])/10
        if l.split(',')[7][-1] == 'E':
            lon = lon
        else:
            lon = -lon
        latt.append(lat)
        lont.append(lon)
        time.append(str(int(l.split(',')[2])))
        intt.append(float(l.split(',')[8]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    intt = np.asarray(intt)
    time_track, ind = np.unique(time,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]
    int_track = intt[ind]
    name_storm = l.split(',')[-11]

    return lon_track, lat_track, time_track, int_track, name_storm

################################################################################
def get_GFS_track_and_int(file_GFS_track,cycle):

    ff = open(file_GFS_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    lead_time = []
    cycl = []
    intt = []
    for l in f:
        if np.logical_and(str(int(l.split(',')[2])) == cycle, \
                      l.split(',')[4].split()[0] == 'AVNO'):

            lat = float(l.split(',')[6][0:4])/10
            if l.split(',')[6][-1] == 'N':
                lat = lat
            else:
                lat = -lat
            lon = float(l.split(',')[7][0:5])/10
            if l.split(',')[7][-1] == 'E':
                lon = lon
            else:
                lon = -lon
            latt.append(lat)
            lont.append(lon)
            lead_time.append(int(l.split(',')[5][1:4]))
            cycl.append(str(int(l.split(',')[2])))
            intt.append(float(l.split(',')[8]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    lead_time_track = np.asarray(lead_time)
    cycle_t = np.asarray(cycl)
    intt = np.asarray(intt)
    lead_time_track, ind = np.unique(lead_time_track,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]
    int_track = intt[ind]
    cycle_track = cycle_t[ind]

    return lon_track, lat_track, lead_time_track, int_track, cycle_track

################################################################################
#%% Conversion from geographic longitude and latitude to HYCOM convention
def geo_coord_to_HYCOM_coord(long,latg):

    long = np.asarray(long)
    if np.ndim(long) > 0:
        lonm = [360 + ln if ln<0 else ln for ln in long]
    else:
        lonm = [360 + long if long<0 else long][0]
    latm = latg

    return lonm, latm

################################################################################
#%% Conversion from geographic longitude and latitude to HYCOM convention
def HYCOM_coord_to_geo_coord(lonh,lath):
    lonh = np.asarray(lonh)
    if np.ndim(lonh) > 0:
        long = [ln-360 if ln>=180 else ln for ln in lonh]
    else:
        long = [lonh-360 if lonh>=180 else lonh][0]
    latg = lath
    return long, latg

################################################################################
def geo_coor_to_GOFS_coord(long,latg):

    target_lon = np.empty((len(long),))
    target_lon[:] = np.nan
    for i,ii in enumerate(long):
        if ii < 0:
            target_lon[i] = 360 + ii
        else:
            target_lon[i] = ii
    target_lat = latg

    return target_lon, target_lat

################################################################################
def GOFS_coor_to_geo_coord(lon_GOFS,lat_GOFS):

    lon_geo = np.empty((len(lon_GOFS),))
    lon_geo[:] = np.nan
    for i in range(len(lon_GOFS)):
        if lon_GOFS[i] > 180:
            lon_geo[i] = lon_GOFS[i] - 360
        else:
            lon_geo[i] = lon_GOFS[i]
    lat_geo = lat_GOFS

    return lon_geo, lat_geo

################################################################################
def Haversine(lat1,lon1,lat2,lon2):
    """
    This uses the haversine formula to calculate the great-circle distance
    between two points that is, the shortest distance over the earth surface
    giving an the-crow-flies~@~Y distance between the points
    (ignoring any hills they fly over, of course!).
    Haversine formula:
    a = sin²(~T~F/2) + cos ~F1 ~K~E cos ~F2 ~K~E sin²(~Tλ/2)
    c = 2 ~K~E atan2( ~H~Za, ~H~Z(1~H~Ra) )
    d = R ~K~E c
    where   ~F is latitude, λ is longitude, R is earth's radius
    (mean radius = 6,371km). Note that angles need to be in radians to pass to
    trig functions!
    """

    R = 6371.0088
    lat1,lon1,lat2,lon2 = map(np.radians, [lat1,lon1,lat2,lon2])

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2) **2
    c = 2 * np.arctan2(a**0.5, (1-a)**0.5)
    d = R * c

    tan_theta = dlat/dlon
    cos_theta = (dlon*180/np.pi)*111.111/d
    sin_theta = (dlat*180/np.pi)*111.111/d

    return d, tan_theta, cos_theta, sin_theta

################################################################################
def get_glider_transect_from_HAFS_HYCOM(ncfiles,lon,lat,depth,tstamp_glider,long,latg):

    target_temp = np.empty((len(depth),len(ncfiles)))
    target_temp[:] = np.nan
    target_salt = np.empty((len(depth),len(ncfiles)))
    target_salt[:] = np.nan
    target_depth = np.empty((len(depth),len(ncfiles)))
    target_depth[:] = np.nan
    time = []

    for x,file in enumerate(ncfiles):
        print(x)
        model = xr.open_dataset(file)
        t = model['MT'][:]
        timestamp = mdates.date2num(t)[0]
        time.append(mdates.num2date(timestamp))

        # Interpolating latg and longlider into HYCOM grid
        sublon = np.interp(timestamp,tstamp_glider,long)
        sublat = np.interp(timestamp,tstamp_glider,latg)
        oklon = np.int(np.round(np.interp(sublon,lon,np.arange(len(lon)))))
        oklat = np.int(np.round(np.interp(sublat,lat,np.arange(len(lat)))))

        target_temp[:,x] = np.asarray(model['temperature'][0,:,oklat,oklon])
        target_salt[:,x] = np.asarray(model['salinity'][0,:,oklat,oklon])

    return time, target_temp, target_salt

################################################################################
def get_var_from_model_following_trajectory(files_model,var_name,time_name,lon_name,lat_name,lon_obs,lat_obs,timestamp_obs,**kwargs):

    # depth_level is an optional argument
    # kwargs = dict(depth_level = 0)
 
    depth_level = kwargs.get('depth_level', None)

    target_var = np.empty((len(files_model)))
    target_var[:] = np.nan
    target_time = []

    for x,file in enumerate(files_model):
        print(x)
        model = xr.open_dataset(file,engine="pynio")
        if file.split('.')[-1] == 'nc':
            t = model[time_name][:]
            timestamp = mdates.date2num(t)[0]
            target_time.append(mdates.num2date(timestamp))
            lon_model = np.asarray(model[lon_name][:])
            lat_model = np.asarray(model[lat_name][:])
        if file.split('.')[-1] == 'grb2':
            t0 = model['TMP_P0_L1_GLL0'].attrs['initial_time']
            dt = model['TMP_P0_L1_GLL0'].attrs['forecast_time'][0]
            t = datetime.strptime(t0, '%m/%d/%Y (%H:%M)') + timedelta(hours=int(dt))
            timestamp = mdates.date2num(t)
            target_time.append(t)
            lon_model = np.asarray(model.lon_0)
            lat_model = np.asarray(model.lat_0)
        
        # Interpolating latg and longlider into HYCOM grid
        sublon = np.interp(timestamp,timestamp_obs,lon_obs)
        sublat = np.interp(timestamp,timestamp_obs,lat_obs)
        oklon = int(np.round(np.interp(sublon,lon_model,np.arange(len(lon_model)))))
        oklat = int(np.round(np.interp(sublat,lat_model,np.arange(len(lat_model)))))
        if file.split('.')[-1] == 'nc':
            if model[var_name].ndim == 4:
                target_var[x] = np.asarray(model[var_name][0,depth_level,oklat,oklon])
            if model[var_name].ndim == 3:
                target_var[x] = np.asarray(model[var_name][0,oklat,oklon])
        if file.split('.')[-1] == 'grb2':
            if model[var_name].ndim == 3:
                target_var[x] = np.asarray(model[var_name][depth_level,oklat,oklon])
            if model[var_name].ndim == 2:
                target_var[x] = np.asarray(model[var_name][oklat,oklon])

    return target_time, target_var

################################################################################
def figure_transect_temp(time,depth,temp,date_ini,date_end,max_depth,kw,color_map):

    #Time window
    year_ini = int(date_ini.split('/')[0])
    month_ini = int(date_ini.split('/')[1])
    day_ini = int(date_ini.split('/')[2])

    year_end = int(date_end.split('/')[0])
    month_end = int(date_end.split('/')[1])
    day_end = int(date_end.split('/')[2])

    tini = datetime(year_ini, month_ini, day_ini)
    tend = datetime(year_end, month_end, day_end)

    fig, ax = plt.subplots(figsize=(7, 3))
    #cs = plt.contourf(time,depth,temp,cmap=cmocean.cm.thermal,**kw)
    cs = plt.contourf(time,depth,temp,cmap=color_map,**kw)
    plt.contour(time,depth,temp,levels=[26],colors = 'k')
    ax.set_xlim(time[0], time[-1])
    ax.set_ylabel('Depth (m)',fontsize=14)
    cbar = plt.colorbar(cs)
    cbar.ax.set_ylabel('($^oC$)',fontsize=14)

    xvec = [tini + timedelta(int(dt)) for dt in np.arange((tend-tini).days+1)[::2]]
    plt.xticks(xvec,fontsize=12)
    xfmt = mdates.DateFormatter('%b-%d')
    ax.xaxis.set_major_formatter(xfmt)
    plt.ylim(-np.abs(max_depth),0)
    plt.xlim(tini,tend)

################################################################################
def glider_data_vector_to_array(depth,time,var,lat,lon):

    upcast = np.where(np.diff(depth) < 0)[0]
    oku = np.where(np.diff(upcast)>1)[0]
    end_upcast = upcast[oku]

    downcast = np.where(np.diff(depth) > 0)[0]
    okd = np.where(np.diff(downcast)>1)[0]
    end_downcast = downcast[okd]

    inn = np.hstack([0,np.unique(np.hstack([end_upcast,end_downcast])),len(depth)])
    okp = np.where(np.diff(inn) > 10)[0]
    ind = np.hstack([inn[0],inn[okp+1]])
    zn = np.max(np.diff(ind)) + 1

    depthg = np.empty((zn,len(ind)))
    depthg[:] = np.nan
    timeg = np.empty((zn,len(ind)))
    timeg[:] = np.nan
    varg = np.empty((zn,len(ind)))
    varg[:] = np.nan
    latg = np.empty((zn,len(ind)))
    latg[:] = np.nan
    long = np.empty((zn,len(ind)))
    long[:] = np.nan

    for i in np.arange(len(ind)):
        if i == 0:
            #print(i)
            indd = np.argsort(depth[ind[i]:ind[i+1]+2])
            depthg[0:len(depth[ind[i]:ind[i+1]+2]),i] = depth[ind[i]:ind[i+1]+2][indd]
            timeg[0:len(depth[ind[i]:ind[i+1]+2]),i] = mdates.date2num(time[ind[i]:ind[i+1]+2][indd])
            varg[0:len(var[ind[i]:ind[i+1]+2]),i] = var[ind[i]:ind[i+1]+2][indd]
            latg[0:len(var[ind[i]:ind[i+1]+2]),i] = lat[ind[i]:ind[i+1]+2][indd]
            long[0:len(var[ind[i]:ind[i+1]+2]),i] = lon[ind[i]:ind[i+1]+2][indd]
        if i < len(ind)-1:
            #print(i)
            indd = np.argsort(depth[ind[i]+1:ind[i+1]+2])
            depthg[0:len(depth[ind[i]+1:ind[i+1]+2]),i] = depth[ind[i]+1:ind[i+1]+2][indd]
            timeg[0:len(depth[ind[i]+1:ind[i+1]+2]),i] = mdates.date2num(time[ind[i]+1:ind[i+1]+2][indd])
            varg[0:len(var[ind[i]+1:ind[i+1]+2]),i] = var[ind[i]+1:ind[i+1]+2][indd]
            latg[0:len(var[ind[i]+1:ind[i+1]+2]),i] = lat[ind[i]+1:ind[i+1]+2][indd]
            long[0:len(var[ind[i]+1:ind[i+1]+2]),i] = lon[ind[i]+1:ind[i+1]+2][indd]
        else:
            #print(i)
            indd = np.argsort(depth[ind[i]+1:len(depth)])
            depthg[0:len(depth[ind[i]+1:len(depth)]),i] = depth[ind[i]+1:len(depth)][indd]
            timeg[0:len(depth[ind[i]+1:len(depth)]),i] = mdates.date2num(time[ind[i]+1:len(depth)][indd])
            varg[0:len(var[ind[i]+1:len(var)]),i] = var[ind[i]+1:len(var)][indd]
            latg[0:len(var[ind[i]+1:len(var)]),i] = lat[ind[i]+1:len(var)][indd]
            long[0:len(var[ind[i]+1:len(var)]),i] = lon[ind[i]+1:len(var)][indd]

    return depthg, timeg, varg, latg, long

################################################################################
def grid_glider_data(varg,timeg,depthg,delta_z=0.3):

    # sort time variable
    okt = np.argsort(timeg)
    timeg_gridded = timeg[okt]
    depthgg = depthg[:,okt]
    vargg = varg[:,okt]

    # Grid variables
    depthg_gridded = np.arange(0,np.nanmax(depthgg),delta_z)
    varg_gridded = np.empty((len(depthg_gridded),len(timeg)))
    varg_gridded[:] = np.nan

    for t,tt in enumerate(timeg):
        depthu,oku = np.unique(depthgg[:,t],return_index=True)
        varu = vargg[oku,t]
        okdd = np.isfinite(depthu)
        depthf = depthu[okdd]
        varf = varu[okdd]
        ok = np.asarray(np.isfinite(varf))
        if np.sum(ok) < 3:
            print('sum less 3 ',t)
            varg_gridded[:,t] = np.nan
        else:
            #okd = depthg_gridded < np.max(depthf[ok])
            okd = np.logical_and(depthg_gridded >= np.min(depthf[ok]),\
                                 depthg_gridded < np.max(depthf[ok]))
            varg_gridded[okd,t] = np.interp(depthg_gridded[okd],depthf[ok],varf[ok])

    return varg_gridded, timeg_gridded, depthg_gridded
################################################################################

