
#%% User input

scratch_dir = '/scratch2/NOS/nosofs/Maria.Aristizabal/'

#model = 'HWRF2019_POM'
#mtype = 'oper'

#model = 'HWRF2020_POM'
#mtype = 'exp'

model = 'HWRF2020_HYCOM'
mtype = 'exp'

storm_name = 'Dorian'
storm_id = 'dorian05l'
#cycles = ['2019082800','2019082806','2019082812','2019082818','2019082900','2019082906','2019082912','2019082918','2019083000','2019083006','2019083012','2019083018','2019083100','2019083106','2019083112','2019083118','2019090100','2019090106']
cycles = ['2019083012']

import os
import glob

for cycle in cycles:
    print(cycle)
    Dir_target = scratch_dir + model + '_' + storm_name + '/' + model + '_' + storm_id + '.' + cycle + '_' + mtype + '/'
 
    os.chdir(Dir_target)
    os.system('tar -xf ' + storm_id + '.' + cycle + '.tar' + ' ' + storm_id + '.' + cycle + '.' + 'trak.hwrf.atcfunix') 
    #os.system('rm ' + storm_id + '.' + cycle + '.tar')

os.chdir('/home/Maria.Aristizabal/Dorian_2019/Code')
'''
os.system('wgrib2 -v /scratch2/NCEPDEV/stmp1/Hyun.Sook.Kim/HWRF19/com/2019083018/05L/dorian05l.2019083018.hwrfprs.storm.0p015.f003.grb2 -match ":(SHTFL|LHTFL|DSWRF|DLWRF|USWRF|ULWRF|WTMP):(surface)|:(UGRD|VGRD):(10 m above ground)" -netcdf /home/Maria.Aristizabal/Dorian_2019/HWRF_nc_files/dorian05l.2019083018.hwrfprs.storm.0p015.f003.nc')
'''
