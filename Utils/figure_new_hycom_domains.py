
### User input

# hycom_nhc domain: 
lon_nhc = [-178.000, 14.960]
lat_nhc = [-23.033, 45.494]

# hycom_jtnh domain: 
lon_jtnh = [30, 179.920]
lat_jtnh = [-23.033, 45.494]

# hycom_jtnh domain: 
lon_jtsh = [24.960, 120.000]
lat_jtsh = [-45.987, 7.022]

domains = ['hycom_nhc','hycom_jtnh','hycom_jtsh']

bath_file = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/bathymetry_files/GEBCO_2021_sub_ice_topo.nc'

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

def coast180():
   inf=open('/home/Maria.Aristizabal/fixdata/coastlines_180x180.dat','r')
   c1 = []
   c2 = []
   hsk=np.genfromtxt(inf)
   c1=hsk[:,0]
   c2=hsk[:,1]
   return (c1,c2)

#################################################################################
##% Read coast line
cx,cy = coast180()

#################################################################################
#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = np.asarray(ncbath.variables['lat'][::4])
bath_lon = np.asarray(ncbath.variables['lon'][::4])
bath_elev = np.asarray(ncbath.variables['elevation'][::4,::4])

#################################################################################

kw = dict(levels=np.arange(-8000,8001,100)) 
fig,ax = plt.subplots(figsize=(9,5))
plt.contourf(bath_lon,bath_lat,bath_elev,cmap='GnBu_r',**kw)
plt.colorbar()
plt.plot(cx,cy,'-',color='gray',markersize=2)
plt.axis('scaled')

for domain in domains:
    if domain == 'hycom_nhc':
        lon = lon_nhc
        lat = lat_nhc
    if domain == 'hycom_jtnh':
        lon = lon_jtnh
        lat = lat_jtnh
    if domain == 'hycom_jtsh':
        lon = lon_jtsh
        lat = lat_jtsh
    
    lon_vec = np.arange(lon[0],lon[1])
    lat_vec = np.arange(lat[0],lat[1])
    
    plt.plot(lon_vec,np.tile(lat[0],len(lon_vec)),'-k')
    plt.plot(lon_vec,np.tile(lat[1],len(lon_vec)),'-k')
    plt.plot(np.tile(lon[0],len(lat_vec)),lat_vec,'-k')
    plt.plot(np.tile(lon[1],len(lat_vec)),lat_vec,'-k')

plt.ion()
plt.show()
plt.xlim([-180,180])

plt.savefig('hycom_new_domains2',dpi=600)
