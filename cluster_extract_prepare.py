import numpy as np
import nibabel as nb
import seaborn as sb
import pandas as pd
import glob
from scipy.ndimage import measurements
import os
import copy


'''
Reads in dataframe output from prepare_clinical.py, add data on if participant has lesion mask overlap from VLSM mask associated with increased MADRS score, extract FALFF values from clusters found to show an interaction and save for analysis in R.
'''

subj_dir = '/home/peter/Desktop/prepare/rest/output/'

#Load spreadsheet output from analysis script
df = pd.read_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs.csv', sep = ',')

##Lesion mask overlap with VLSM results.##


#df['mask_vlsm_overlap'] = np.zeros(len(df))

#Collect lesion mask files.
#lesion_mask_dir = subj_dir + 'lesion_masks/'
#lesion_masks = os.listdir(lesion_mask_dir + 'baseline/')
#lesion_masks.sort()

#Load in significant voxels from VLSM analysis.
#vlsm_mask = nb.load(lesion_mask_dir + 'madrs_vlsm_thresh.nii').get_data()
#vlsm_mask[np.isnan(vlsm_mask) == 1] = 0 #Zero out nans.


#Loop - read in mask, check the initials in the mask filename are in the corresponding row of the data sheet.
#Check if lesion mask overlaps with VLSM sig voxels and add a 1 to mask_vlsm_overlap if they do, else 0 if not.
#for n, mask in enumerate(lesion_masks):
#    if df['initials'].iloc[n] in mask:
#        mask_data = nb.load(lesion_mask_dir + 'baseline/' + mask).get_data()
#        mask_data[np.isnan(mask_data) == 1] = 0 #Zero out nans.
#        overlap_mask = vlsm_mask + mask_data
#        if overlap_mask.max() == 2:
#            #df['mask_vlsm_overlap'].iloc[n] = 1
#            df['mask_vlsm_overlap'].iloc[n] = sum((vlsm_mask + mask_data) == 2)
#        else:
#             df['mask_vlsm_overlap'].iloc[n] = 0
#    else:
#        print('\nWARNING! THE INITIALS {} ARE NOT IN {}! Exiting.'.format(df['initials'].iloc[n], mask))
#        break

#Save output
#df.to_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs.csv', sep = ',', index = False)        
    
    



##fALFF data extraction.##

thresh_dir = '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/'


#Non range standardised.
cond_dir_list = [thresh_dir + 'falff_analysis_noglob', thresh_dir + 'falff_slow_4_analysis_noglob', thresh_dir + 'falff_slow_5_analysis_noglob']

#Range standardised.
#cond_dir_list = [thresh_dir + 'falff_analysis_stand', thresh_dir + 'falff_slow_4_analysis_stand']

cond_fn = ['szfalff_av-noglobal_mni.nii', 'szfalff_slow_4-noglobal_mni.nii', 'szfalff_slow_5-noglobal_mni.nii' ]

for n, cond_dir in enumerate(cond_dir_list):
    
    cond = os.path.split(cond_dir)[1]
    sig_dir = cond_dir + '/l2analysis/EstimateContrasts/'
    label_output_dir = sig_dir + '/labelled_clusters/'
 
    #Non range standardised.
    
    
    #Get paricipant files to extract data.
    df[cond + '_fn'] = np.zeros(len(df)) #Create padded column
    falff_files = glob.glob(subj_dir + '/*/f_alff/' + cond_fn[n])
    falff_files =  [x for x in falff_files if '12months' not in x] #Drop 12 month files
    falff_files.sort()
    for n, fn in enumerate(falff_files):
        if df['initials'].iloc[n] in fn:
            df[cond + '_fn'].iloc[n] = fn
        else:
            print('\nWARNING! THE INITIALS {} ARE NOT IN {}! Exiting.'.format(df['initials'].iloc[n], fn))
            break

#SPM uses an 18 sided structural element to determine individual clusters. 
#18 connections
    strct = np.ones((3, 3, 3))
    strct[0, 0, 0] = 0
    strct[0, 2, 0] = 0
    strct[2, 0, 0] = 0
    strct[2, 2, 0] = 0
    strct[0, 0, 2] = 0
    strct[0, 2, 2] = 0
    strct[2, 0, 2] = 0
    strct[2, 2, 2] = 0
     
        
    try:
        os.mkdir(label_output_dir)
    except:
        print('Output directory exists. Check there first...')
    
    
    #load_file = glob.glob(sig_dir + '/falff*.nii')[0]
    load_file = glob.glob(sig_dir + '/prepare_hGTl_thresh.nii')[0]     
    label_info = nb.load(load_file)
    label_data = label_info.get_data()
    label_data[np.isnan(label_data) == 1] = 0 #Un-nan data
    label_data[label_data > 0] = 1 #Binarise
    
    #Label individual clusters in thresh_data using 18 sided element, save nii to compare output with threshold map.
    lw, num = measurements.label(label_data,structure = strct)
    
    img_filename = label_output_dir + os.path.split(load_file)[1][:-4] + '_clust_label.nii'
    img = nb.Nifti1Image(lw, header = label_info.header, affine = label_info.affine)
    img.to_filename(img_filename)
    
    
    #Loop over unique values in labelled threshold mask
    for clustn in np.arange(0,num):
        idx = clustn + 1
        mask = lw == idx #Make a binary mask of only the cluster values currently on loop.
        
        df[cond + '_clustval_' + str(idx)] = np.zeros(len(df))
        print('\nExtracting falff values for {} condition in cluster {}'.format(cond, idx))
        
    
    #Read subjected nii in, make 1d array of data
        for i, niifile in enumerate(df[cond + '_fn']):
            print('\nStarting extraction on subject {}'.format((niifile,)))
            nii = nb.load(niifile).get_data()            
            clust_data = nii[mask]
    #Collect voxel value from labelled sig region, write to dataframe     
            clustmean = np.mean(nii[mask])
            df[cond + '_clustval_'+str(idx)][i] = clustmean
    
#Save to csvstroke
#if df['mask_vlsm_overlap'].max() == 1:
#    df.to_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_vlsm_falff.csv')
#else:
#    df.to_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_falff.csv')
df.to_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_falff.csv')   


#Plot extracted data
plot_df = copy.deepcopy(df)

#Rename variables
plot_df.columns.values[plot_df.columns == 'MADRS_score_3mth'] = 'MADRS Score'
plot_df.columns.values[plot_df.columns == 'madrs_group'] = 'Group'
plot_df.replace(to_replace = {'Group': {0 : 'Low', 1 : 'High'}}, inplace = True)

clust_list = ['falff_analysis_noglob_clustval_1', 
'falff_analysis_noglob_clustval_2', 
'falff_slow_4_analysis_noglob_clustval_1', 
'falff_slow_4_analysis_noglob_clustval_2', 
'falff_slow_5_analysis_noglob_clustval_1']

region_names = ['fALFF Broadband Left Superior Temporal', 
'fALFF Broadband Left Insula', 
'fALFF Slow 4 Left Thalamus', 
'fALFF Slow 4 Right Caudate', 
'fALFF Slow 5 Left Cerebellum (Posterior Lobe)']

var_list = zip(clust_list, region_names)
markers = ['o', 'x']

#Set figure properties
sb.set_context('talk', font_scale = 2)
sb.set_style('dark')

#Plot figures
for var in var_list:
    plot_df.rename(columns = {var[0] : var[1]}, inplace = True)
    #fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    sb.lmplot(x = var[1], y = 'MADRS Score', hue = 'Group', data = plot_df, markers = markers, scatter_kws = {'s' : 50}, line_kws = {'alpha' : 0.50, 'lw' : 1}, palette = {'Low' : 'g', 'High' : 'b'}, legend = False, aspect = 2)
    sb.plt.ylim([-15, 40])
    plt.legend(loc = 'upper left')
    image_fn = '/home/peter/Dropbox/post_doc/florey_leeanne/papers/prepare/madrs_baseline/images/' + var[1].replace(' ', '_').lower() + '_interaction'
    plt.savefig(image_fn, format = 'png')
    plt.close()


#g = [sb.lmplot(x = var, y = 'MADRS_score_3mth', hue = 'madrs_group', data = df) for var in var_list]
#g2 = [sb.lmplot(x = var, y = 'madrs_score_stand', hue = 'madrs_group', data = df) for var in var_list]


