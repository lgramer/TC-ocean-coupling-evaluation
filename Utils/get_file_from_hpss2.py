
#%% User input

# DO THIS IN THE TERMINAL BEFORE RUNNING THIS SCRIPT!!!!!!!!
#%% Load module to access the hpss system
#os.system('module load hpss')
#%% Load neccesary modules to use grib2 
#os.system('module load intel/19.0.5.281 wgrib2/2.0.8')

#model = 'HWRF2019_POM'
#Dir_hpss = '/NCEPPROD/2year/hpssprod/runhistory/rh2019/hwrf/05l/'
#mtype = 'oper'

#model = 'HWRF2020_POM'
#mtype = 'exp'
#Dir_hpss = '/NCEPDEV/emc-hwrf/2year/emc.hurpara/H20C/'

model = 'HWRF2020_HYCOM'
mtype = 'exp'
Dir_hpss = '/NCEPDEV/emc-hwrf/2year/emc.hurpara/H20H/'

storm_name = 'Dorian'
storm_id = 'dorian05l'
#cycles = ['2019082806','2019082812','2019082818','2019082900','2019082906','2019082912','2019082918','2019083000','2019083006','2019083012','2019083018','2019083100','2019083106','2019083112','2019083118','2019090100','2019090106']
#cycles = ['2019083106','2019083112','2019083118','2019090100','2019090106']
cycles = ['2019083012']

import os
import glob

for cycle in cycles:
    print(cycle)
    Dir_target = '/scratch2/NOS/nosofs/Maria.Aristizabal/' + model + '_' + storm_name + '/' + model + '_' + storm_id + '.' + cycle + '_' + mtype + '/'
    Dir_output = '/scratch2/NOS/nosofs/Maria.Aristizabal/' + model + '_' + storm_name + '/' +model + '_' + storm_id + '.' + cycle + '_grb2_to_nc_' + mtype + '/'
 
    os.system('mkdir ' + Dir_target)
    os.system('mkdir ' + Dir_output)
    os.chdir(Dir_target)
    os.system('hsi get ' + Dir_hpss + storm_id + '.' + cycle + '.tar')
    os.system('tar -xf ' + storm_id + '.' + cycle + '.tar' + ' *hwrfprs*storm*grb2*') 

    #%% Find HWRF files
    HWRF_files = sorted(glob.glob(os.path.join(Dir_target,'*hwrfprs*storm*grb2')))

    #%% convert from  grib2 to nc files
    for fl in HWRF_files:
        #print(fl)
        ll = 'wgrib2 -v ' + fl + ' -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf ' + Dir_output + fl.split('/')[-1][0:-4] +'nc'

        os.system(ll) 
        
    #os.system('rm ' + storm_id + '.' + cycle + '.tar')

os.chdir('/home/Maria.Aristizabal/Dorian_2019/Code')

'''
os.system('wgrib2 -v /scratch2/NCEPDEV/stmp1/Hyun.Sook.Kim/HWRF19/com/2019083018/05L/dorian05l.2019083018.hwrfprs.storm.0p015.f003.grb2 -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf /home/Maria.Aristizabal/Dorian_2019/HWRF_nc_files/dorian05l.2019083018.hwrfprs.storm.0p015.f003.nc')
'''
