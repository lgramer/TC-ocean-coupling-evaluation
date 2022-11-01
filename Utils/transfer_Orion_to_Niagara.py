import os, sys, subprocess, time

#os.system('module load globus-cli')
#os.system('export ORION=84bad22e-cb80-11ea-9a44-0255d23c44ef')
#os.system('export NIAGARA=21467dd0-afd6-11ea-8f12-0a21f750d19b')
#os.system('globus endpoint activate --web --no-browser $ORION')
#os.system('globus endpoint activate --web --no-browser $NIAGARA')

orion_dir = '$ORION\:/work/noaa/hwrf/noscrub/maristiz/hafsarch/hafsv0p2a_phase3/'
niagara_dir = '$NIAGARA\:/collab1/data/Maria.Aristizabal/hafsv0p2a_phase3/'
hpss_dir = '/5year/NCEPDEV/emc-hwrf/Maria.Aristizabal/hafsv0p2a_phase3/'

niagara_dir2 = niagara_dir.split(':')[1]

globus_ls_cmd = 'globus ls '+ orion_dir + ' --filter \'~*.done\''

result = subprocess.check_output(globus_ls_cmd,shell=True)

#for donefile in result.split('\n'):
for donefile in result.split():
    if str(donefile)[:-1].endswith('done'):
        tarfile = str(donefile)[2:-6]
# transfer from orion to niagara
        globus_transfer_cmd = 'globus transfer --notify off ' + orion_dir + tarfile + ' ' + niagara_dir+tarfile
#        globus_transfer_cmd='globus transfer '+orion_dir+tarfile+' '+niagara_dir+tarfile+' --notify failed,inactive'
#        globus_transfer_cmd="globus transfer "+orion_dir+tarfile+" "+niagara_dir+tarfile+" --notify 'off'"
#        globus_transfer_cmd='globus transfer '+orion_dir+tarfile+' '+niagara_dir+tarfile
        print (globus_transfer_cmd)
        result2 = subprocess.check_output(globus_transfer_cmd,shell=True)
        for line in str(result2).split('\\n'):
            if (line.startswith('Task ID')) :
                task_id = line.split(': ')[-1]
# check globus transfer status
        globus_check_cmd = 'globus task show ' + task_id
        job_status = 1
        for i in range(1,61):
            result3 = subprocess.check_output(globus_check_cmd,shell=True)
            for line in str(result3).split('\\n'):
                if (line.startswith('Status')):
                    task_status = line[-9:]
                    print (task_status)
            if (task_status == 'SUCCEEDED'):
                job_status = 0
                break
            time.sleep(60)
        if (job_status == 0):
# delete orion files
            globus_delete_cmd1 = 'globus delete --notify off ' + orion_dir + tarfile
            globus_delete_cmd2 = 'globus delete --notify off ' + orion_dir + str(donefile)[2:-1]
            result4 = subprocess.check_output(globus_delete_cmd1,shell=True)
            result5 = subprocess.check_output(globus_delete_cmd2,shell=True)
# push to hpss
            hsi_command = 'hsi put ' + niagara_dir2+tarfile +' : ' + hpss_dir+tarfile
            print (hsi_command)
            os.system(hsi_command)
            if os.path.exists(niagara_dir2+tarfile):
                os.remove(niagara_dir2+tarfile)
            else:
                print(niagara_dir + tarfile+" does not exist")

print ('all done!')
