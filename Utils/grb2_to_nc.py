
#%% User input

#cycle ='2020082512'
#cycle ='2020091100'
cycle ='2019082800'

#Dir_input = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.' + cycle + '_oper/'
#Dir_output = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2019_POM_Dorian/HWRF2019_POM_dorian05l.' + cycle + '_grb2_to_nc_oper/'

Dir_input = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.' + cycle + '_exp/'
Dir_output = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_POM_Dorian/HWRF2020_POM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

#Dir_input = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_HYCOM_Dorian/HWRF2020_HYCOM_dorian05l.' + cycle + '_exp/'
#Dir_output = '/scratch2/AOML/aoml-phod/Maria.Aristizabal/HWRF2020_HYCOM_Dorian/HWRF2020_HYCOM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

#Dir_input = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/HAFSv0p2a_phase3_final/'+cycle+'/00L/'
#Dir_output = Dir_input

import os 
import glob

#%% Load neccesary modules to use grib2 
#os.system('module load intel/19.0.5.281 wgrib2/2.0.8')
#%% Find files
grb2_files = sorted(glob.glob(os.path.join(Dir_input,'*hwrfprs*storm*grb2')))
#grb2_files = sorted(glob.glob(os.path.join(Dir_input,'*hafsprs*synoptic*grb2')))


#%% convert from  grib2 to nc files
for fl in grb2_files[::2][1:3]:
	print(fl)
	#ll = 'wgrib2 -v ' + fl + ' -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP|PRES|TMP):(surface)|:(UGRD|VGRD):(10 m above ground)|:(RH|TMP):(2 m above ground)" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'nc'
	ll = 'wgrib2 -v ' + fl + ' -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP|PRES|TMP):(surface)|:(UGRD|VGRD):(10 m above ground)|:(RH|TMP|UGRD|VGRD)" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'nc'

	os.system(ll) 

'''
#%% convert from  grib2 to nc files
for fl in grb2_files:
	print(fl)
	ll = 'wgrib2 -v ' + fl + ' -match ":(UGRD|VGRD):850 mb|:(UGRD|VGRD):200 mb:" -netcdf ' + Dir_output + fl.split('/')[-1][0:-5] + '_shear' + '.nc'

	os.system(ll) 
'''
'''
os.system('wgrib2 -v /scratch2/NCEPDEV/stmp1/Hyun.Sook.Kim/HWRF19/com/2019083018/05L/dorian05l.2019083018.hwrfprs.storm.0p015.f003.grb2 -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf /home/Maria.Aristizabal/Dorian_2019/HWRF_nc_files/dorian05l.2019083018.hwrfprs.storm.0p015.f003.nc')
'''
