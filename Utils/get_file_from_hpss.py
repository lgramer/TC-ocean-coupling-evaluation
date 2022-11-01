
#%% User input

# DO THIS IN THE TERMINAL BEFORE RUNNING THIS SCRIPT!!!!!!!!
#%% Load module to access the hpss system
#os.system('module load hpss')
#%% Load neccesary modules to use grib2 
#os.system('module load intel/19.0.5.281 wgrib2/2.0.8')

#scratch_dir = '/scratch2/NOS/nosofs/Maria.Aristizabal/'
scratch_dir = '/scratch1/NCEPDEV/hwrf/noscrub/Maria.Aristizabal/'

model = 'HWRF2021_oper'
Dir_hpss = '/NCEPPROD/2year/hpssprod/runhistory/rh2021/hwrf/05l/'

storm_name = 'elsa'
storm_id = '05l'
cycles = ['2021070200']
#cycles = ['2021070106','2021070112','2021070118','2021070206','2021070212','2021070218','2021070300']

import os
import glob

for cycle in cycles:
    print(cycle)
    Dir_target = scratch_dir + model + '/' + storm_id + '/' + cycle + '/'
 
    os.system('mkdir ' + Dir_target)
    os.chdir(Dir_target)
    os.system('hsi get ' + Dir_hpss + storm_name + storm_id + '.' + cycle + '.tar')
    os.system('tar -xf ' + storm_name + storm_id + '.' + cycle + '.tar -v --wildcards' + ' *hwrfprs*storm*grb2*') 
    os.system('tar -xf ' + storm_name + storm_id + '.' + cycle + '.tar -v --wildcards' + ' *pom*00*nc*') 
    os.system('tar -xf ' + storm_name + storm_id + '.' + cycle + '.tar -v --wildcards' + ' *pom*grid*nc*') 
    os.system('tar -xf ' + storm_name + storm_id + '.' + cycle + '.tar -v --wildcards' + ' *trak*hwrf*atcfunix') 

    #%% Find HWRF files
    HWRF_files = sorted(glob.glob(os.path.join(Dir_target,'*hwrfprs*storm*grb2')))

    #%% convert from  grib2 to nc files
    for fl in HWRF_files:
        #print(fl)
        ll = 'wgrib2 -v ' + fl + ' -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf ' + Dir_target + fl.split('/')[-1][0:-4] +'nc'

        os.system(ll) 

    #os.system('rm ' + storm_name + storm_id + '.' + cycle + '.tar')

os.chdir('/home/Maria.Aristizabal/Utils/')

'''
os.system('wgrib2 -v /scratch2/NCEPDEV/stmp1/Hyun.Sook.Kim/HWRF19/com/2019083018/05L/dorian05l.2019083018.hwrfprs.storm.0p015.f003.grb2 -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf /home/Maria.Aristizabal/Dorian_2019/HWRF_nc_files/dorian05l.2019083018.hwrfprs.storm.0p015.f003.nc')
'''
