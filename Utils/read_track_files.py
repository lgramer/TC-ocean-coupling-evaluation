
#%% User input
folder_track_oper = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.2019082800_oper/'
folder_track_exp = '/scratch2/NOS/nosofs/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.2019082800_exp/'
prefix = 'dorian05l.2019082800.trak.hwrf.atcfunix'

#%% 
import numpy as np
from matplotlib import pyplot as plt

#%% Read track files operational
file_track_oper = folder_track_oper + prefix
ff_oper = open(file_track_oper,'r')
f_oper = ff_oper.readlines()

latt_oper = []
lonn_oper = []
intt_oper = []
lead_time_oper = []
for l in f_oper:
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
    latt_oper.append(lat)
    lonn_oper.append(lon)
    intt_oper.append(float(l.split(',')[8]))
    lead_time_oper.append(int(l.split(',')[5][1:4]))

latt_oper = np.asarray(latt_oper)
lonn_oper = np.asarray(lonn_oper)
intt_oper = np.asarray(intt_oper)
lead_time_track_oper, ind = np.unique(lead_time_oper,return_index=True)
lat_track_oper = latt_oper[ind]
lon_track_oper = lonn_oper[ind]
int_track_oper = intt_oper[ind]

#%% Read track files experimental
file_track_exp = folder_track_exp + prefix
ff_exp = open(file_track_exp,'r')
f_exp = ff_exp.readlines()

latt_exp = []
lonn_exp = []
lead_time_exp = []
for l in f_exp:
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
    latt_exp.append(lat)
    lonn_exp.append(lon)
    lead_time_exp.append(int(l.split(',')[5][1:4]))

latt_exp = np.asarray(latt_exp)
lonn_exp = np.asarray(lonn_exp)
lead_time_track_exp, ind = np.unique(lead_time_exp,return_index=True)
lat_track_exp = latt_oper[ind]
lon_track_exp = lonn_oper[ind]


plt.figure()
plt.ion()
plt.plot(lon_track_oper,lat_track_oper,'o-',color='mediumorchid')
plt.plot(lon_track_exp,lat_track_exp,'o-',color='teal')
