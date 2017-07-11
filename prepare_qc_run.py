from __future__ import division
import numpy as np
import nibabel as nb
import glob
import matplotlib.pylab as plt
import pandas as pd
import os
from shutil import copyfile
from shutil import rmtree
import scipy.signal as sig

plt.ioff()

hdir='/home/peter/Desktop/prepare/rest/output/'

subfolds=glob.glob(hdir + 'D_*')
subfolds.sort()

#Make dataframe for housing data.
df = pd.DataFrame(columns=['subs','fd_mean', 'fd_std', 'dvars_mean', 'dvars_std', 'badvols'])

for j,sub in enumerate(subfolds):
    print(os.path.split(sub)[1])

    rotvars = sub + '/QC/rp_aepi_despike.txt'
    
    cleanepi = sub + '/preproc_epis/sbp_residuals_global_trans.nii'
    rawepi = sub + '/preproc_epis/sraepi_despike_trans.nii'
    mmaskfile = sub + '/mni_warped/_apply2mmask0/mmask_trans.nii'
    
    print('Calculating framewise displacement')

    #Calculate framewise displacement:
    motpars = np.loadtxt(rotvars)
    motparsZ = (motpars-np.mean(motpars,axis=0))/np.std(motpars,axis=0)
    
    #compute absolute displacement
    dmotpars = np.zeros(motpars.shape)
    
    dmotpars[1:,:] = np.abs(motpars[1:,:] - motpars[:-1,:])
    
    headradius = 50
    disp = dmotpars.copy()
    disp[:,0:3] = np.pi*headradius*2*(disp[:,0:3]/360)
    
    fd = np.sum(disp,1)
    fd_mean = np.mean(fd)
    fd_std = np.std(fd)


#Calcuate DVARS
    print('Calculating DVARS\n')
    
    epi_clean = nb.load(cleanepi).get_data() #Needed for header and affine info
    epi_raw = nb.load(rawepi).get_data() #Needed for header and affine info
    
    global_mask = nb.load(mmaskfile).get_data()
    global_mask[np.isnan(global_mask)] = 0
    global_mask = global_mask > 0
    

    #Cleaned data
    global_ts_clean = epi_clean[global_mask]
    global_ts_clean_mean = np.mean(global_ts_clean, axis = 0)

    dvars_clean = np.zeros_like(global_ts_clean_mean)
    dvars_clean[1:] = np.sqrt((global_ts_clean_mean[1:] - global_ts_clean_mean[:-1])**2) #Root mean squared difference of volume N to volume N + 1.
    

    #Raw data
    global_ts_raw = epi_raw[global_mask]
    global_ts_raw_mean = np.mean(global_ts_raw, axis = 0)
    dvars_raw = np.zeros_like(global_ts_raw_mean)
    dvars_raw[1:] = np.sqrt((global_ts_raw_mean[1:] - global_ts_raw_mean[:-1])**2) #Root mean squared difference of volume N to volume N + 1.
    
    dvars_clean_mean = np.mean(dvars_clean)
    dvars_clean_std = np.std(dvars_clean)

    badvols=np.where(fd>.5)
    

#Save sub name, framewise displacement and dvars summary stats to dataframe
    df.loc[j,'subs'] = os.path.split(sub)[1]
    df.loc[j,'fd_mean'] = fd_mean
    df.loc[j,'fd_std'] = fd_std
    df.loc[j,'dvars_mean'] = dvars_clean_mean
    df.loc[j,'dvars_std'] = dvars_clean_std
    df.loc[j,'badvols'] = np.shape(badvols)[1]
    
#Horrible, dense, difficult to read patch that plots the data.

    plt.figure(figsize=[25,15],dpi=100)
    plt.subplot(6,1,1); motpars_xyz_plot=plt.plot(motpars[:,:3]);plt.title('Rotation',{'fontsize':10}); plt.ylabel('mm',{'fontsize':10});plt.legend(['x','y','z'],fontsize='xx-small')
    plt.subplot(6,1,2); motpars_rot_plot=plt.plot(motparsZ[:,3:]);plt.title('Translation',{'fontsize':10}); plt.ylabel('rads',{'fontsize':10}); plt.legend(['pitch','yaw','roll'],fontsize='xx-small')
    plt.subplot(6,1,3); fdplot=plt.plot(fd,'r'); plt.plot(np.squeeze(np.zeros([1,len(fd)])+.5,'g')); plt.plot(np.squeeze(np.zeros([1,len(fd)])+.2,'r'));plt.title('Framewise Displacement',{'fontsize':10});plt.ylabel('mm',{'fontsize':10});plt.legend(['FD','.5mm','.2mm'],fontsize='xx-small')
    plt.subplot(6,1,4); dvars_plot=plt.plot(dvars_raw,'r'); plt.title('DVARS raw + clean',{'fontsize':10}); plt.ylabel('Delta Z-score',{'fontsize':10});plt.plot(dvars_clean,'b--');plt.legend(['raw','clean'],fontsize='xx-small')    
    plt.subplot(6,1,5); global_raw_plot=plt.plot(global_ts_raw,'g'); plt.title('Global signal raw',{'fontsize':10}); plt.ylabel('BOLD',{'fontsize':10})
    plt.subplot(6,1,6); global_clean_plot=plt.plot(global_ts_clean,'g'); plt.title('Global signal + clean',{'fontsize':10});plt.ylabel('Zscore'); plt.xlabel('Volumes',{'fontsize':10})
    plt.tight_layout()    
    save_name=os.path.split(sub)[1]+'_qc'
    plt.savefig(sub+'/QC/'+save_name+'.png',format='png')
    print save_name
    plt.close(plt.gcf())

#After all participants are crunched, copy images to output.

try:
    os.mkdir(hdir + 'output_QC')
except:
    print 'Removing output directory'
    rmtree(hdir + 'output_QC')
    
qcfiles=glob.glob(hdir+'*/QC/*_qc.png')
for qclist in qcfiles:
    qcfile=os.path.split(qclist)[1]
    print 'Copying '+qcfile+'to ' + hdir+'output_QC/'
    copyfile(qclist, hdir+'output_QC/'+qcfile)

df.to_csv(hdir + 'prepare_fd_dvars.csv', sep = ',')
