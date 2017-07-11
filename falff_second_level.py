#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Prepare FALFF analysis

Does group level analysis on FALFF data.
Uses healthy controls from the connect dataset and baseline prepare stroke participants.
"""

from nipype.interfaces.spm import OneSampleTTestDesign, TwoSampleTTestDesign, EstimateModel, EstimateContrast
from nipype.pipeline.engine import Workflow, Node
import pandas as pd
import numpy as np
import glob
import os

mask_file = '/home/peter/resources/space_templates/t1_3mm_grey.nii'

maps = ['szalff_av-noglobal_mni.nii', 'szfalff_av-noglobal_mni.nii', 'szfalff_slow_5-noglobal_mni.nii', 'szfalff_slow_4-noglobal_mni.nii'] #Uses z normalised outputs. NOTE: If non-normed data is used, change lines 35 and 37 to reflect change in path length (add 1 digit)

#maps = ['szalff_av-global_mni.nii', 'szfalff_av-global_mni.nii', 'szfalff_slow_5-global_mni.nii', 'szfalff_slow_4-global_mni.nii']

#Read compiled madrs scores (from prepare_clinical output) 
df = pd.read_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs.csv', sep = ',')

df.drop(['Unnamed: 0','Chol_3mth','Trig_3mth', 'HDLC_3mth', 'LDLC_3mth', 'sFol_3mth', 'B12_3mth','CRP_3mth', 'Depress_history_3mths'], axis = 1, inplace = True)

###ONE SAMPLE T (inc. MADRS score).

def prepare_oneT():
    outdirs = ['/home/peter/Desktop/prepare/rest/output/sec_level/one_sample/alff_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/one_sample/falff_analysis_noglob/']
    
    
    for n, outdir in enumerate(outdirs):
    
        stroke_hdir = '/home/peter/Desktop/prepare/rest/output/'
        all_stroke_files = glob.glob(stroke_hdir + '/*/f_alff/' + maps[n])
        stroke_files = [x for x in all_stroke_files if '12months' not in x]
        stroke_files.sort()
        
        #Convert to DataFrame, cull subs that aren't in preprocessed data
        stroke_files = pd.DataFrame(stroke_files, columns = ['files'])
        if n == 0:
            stroke_files['initials'] = os.path.split(os.path.split(stroke_files['files'][n])[0])[0][-2:]
        else:
            stroke_files['initials'] = os.path.split(os.path.split(stroke_files['files'][n])[0])[0][-2:]
        cull_list = stroke_files['initials'].isin(df['initials'])
        df_cull = df[cull_list]
        df_cull.rename(columns = {'Unnamed: 0': 'old_idx'}, inplace = True)
        #Specially remove the second HK (15-28) who refused to die...
        #df_cull = df_cull[df_cull['old_idx'] != 194]
        #df_cull.reset_index(inplace = True, drop = True)
        try:
            np.all(df_cull['initials'] == stroke_files['initials']) == 1       
            stroke_files['madrs_3m'] = df_cull['MADRS_score_3mth']
        except:
            print 'Data are not of identical length. Error in script.'
    
        cov = {'vector': stroke_files['madrs_3m'].tolist(), 'name': 'madrs_3m', 'interaction': 1,  'centering': 1} 
    
        ttest = Node(OneSampleTTestDesign(), name = 'OneSampleT')
        ttest.inputs.in_files = stroke_files['files'].tolist()
        ttest.inputs.covariates = cov
        ttest.inputs.explicit_mask_file = mask_file
        
        modelEst = Node(EstimateModel(), name = 'EstimateModel')
        modelEst.inputs.estimation_method = {'Classical' : 1}
        
        conEst = Node(EstimateContrast(), name = 'EstimateContrasts')
        
        con_1 = ('Stroke','T', ['mean', 'madrs_3m'], [1.0, 0.0])
        con_2 = ('Covariate','T', ['mean','madrs_3m'], [0.0, 1.0])
        
        contrasts = [con_1, con_2]
        
        conEst.inputs.contrasts = contrasts
        conEst.inputs.group_contrast = True
        
        l2analysis=Workflow(name = 'l2analysis')
        l2analysis.base_dir = outdir
            
        l2analysis.connect([(ttest,modelEst,[('spm_mat_file','spm_mat_file')]),
                            (modelEst,conEst,[('spm_mat_file','spm_mat_file'),
                                              ('beta_images','beta_images'),
                                              ('residual_image','residual_image')]),
                                              ])
        
        l2analysis.write_graph(graph2use='colored')
        l2analysis.run('MultiProc', plugin_args={'n_procs':1})



###TWO SAMPLE T###
#Compare low MADRS stroke (<8) with high MADRS stroke (>8) main effects

def prepare_indT_madrs():
    outdirs = ['/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_main/alff_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_main/falff_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_main/falff_slow_5_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_main/falff_slow_4_analysis_noglob/']
    
    
    for n,outdir in enumerate(outdirs):
    
        stroke_hdir = '/home/peter/Desktop/prepare/rest/output/'
        all_stroke_files = glob.glob(stroke_hdir + '/*/f_alff/' + maps[n])
        stroke_files = [x for x in all_stroke_files if '12months' not in x]
        stroke_files.sort()
        
        #Convert to DataFrame, cull subs that aren't in preprocessed data
        stroke_files = pd.DataFrame(stroke_files, columns = ['files'])
        if n == 0:
            stroke_files['initials'] = os.path.split(os.path.split(stroke_files['files'][n])[0])[0][-2:]
        else:
            stroke_files['initials'] = os.path.split(os.path.split(stroke_files['files'][n])[0])[0][-2:]
        cull_list = stroke_files['initials'].isin(df['initials'])
        df_cull = df[cull_list]
        df_cull.rename(columns = {'Unnamed: 0': 'old_idx'}, inplace = True)
        #Specially remove the second HK (15-28) who refused to die...
        #df_cull = df_cull[df_cull['old_idx'] != 194]
        #df_cull.reset_index(inplace = True, drop = True)
        try:
            np.all(df_cull['initials'] == stroke_files['initials']) == 1       
            stroke_files['madrs_3m'] = df_cull['MADRS_score_3mth']
        except:
            print 'Data are not of identical length. Error in script.'
    
    
    #SPLIT GROUPS & MAKE COVARIATE DEFINITIONS
        stroke_low = stroke_files[stroke_files['madrs_3m'] <= 8]
        stroke_high = stroke_files[stroke_files['madrs_3m'] > 8]
    
    
    #INITIALISE T-TEST OBJECT
        ttest = Node(TwoSampleTTestDesign(), name = 'TwoSampleT')
        ttest.inputs.group1_files = stroke_low['files'].tolist()
        ttest.inputs.group2_files = stroke_high['files'].tolist()
        ttest.inputs.explicit_mask_file = mask_file
        
    #Estimate model (aka betas)
        modelEst = Node(EstimateModel(), name = 'EstimateModel')
        modelEst.inputs.estimation_method = {'Classical' : 1}
    
    #Estimate contrasts    
        conEst = Node(EstimateContrast(), name = 'EstimateContrasts')
        
        con_1 = ('mean_low','T', ['Group_{1}', 'Group_{2}'], [1.0, 0.0])
        con_2 = ('mean_high','T', ['Group_{1}','Group_{2}'], [0.0, 1.0])
        con_3 = ('Low>High','T',['Group_{1}','Group_{2}'],[1.0, -1.0])
        con_4 = ('High>Low','T',['Group_{1}','Group_{2}'],[-1.0, 1.0])
    
        contrasts = [con_1, con_2, con_3, con_4]
        
        conEst.inputs.contrasts = contrasts
        conEst.inputs.group_contrast = True
        
    #RUN
        l2analysis=Workflow(name = 'l2analysis')
        l2analysis.base_dir = outdir
            
        l2analysis.connect([(ttest,modelEst,[('spm_mat_file','spm_mat_file')]),
                            (modelEst,conEst,[('spm_mat_file','spm_mat_file'),
                                              ('beta_images','beta_images'),
                                              ('residual_image','residual_image')]),
                                              ])
        
        l2analysis.write_graph(graph2use='colored')
        l2analysis.run('MultiProc', plugin_args={'n_procs':1})






###TWO SAMPLE T###
#Compare low MADRS stroke (<8) with high MADRS stroke (>8) interaction (group x madrs score)
#Note: stand = Standardisation of MADRS scores between 0 and 1
#Note: For modelling the interaction: 
#Interaction = 1 (None), 2 (With Factor 1), 3 (With Factor 2), 4 (With Factor 3)
#centering = 1 (Overall mean), 2 (Factor 1 mean), 3 (Factor 2 mean), 4 (Factor 3 mean), 5 (No Centering), 6 (User Specified Value), 7 (As implied by ANCOVA), 8 (GM)

def prepare_indT_madrs_int(stand):
    outdirs = ['/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/alff_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/falff_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/falff_slow_5_analysis_noglob/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/falff_slow_4_analysis_noglob/']

    #outdirs = ['/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/alff_analysis/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/falff_analysis/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/falff_slow_5_analysis/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_madrs_interaction/falff_slow_4_analysis/']

        
    for n,outdir in enumerate(outdirs):
    
        stroke_hdir = '/home/peter/Desktop/prepare/rest/output/'
        all_stroke_files = glob.glob(stroke_hdir + '/*/f_alff/' + maps[n])
        stroke_files = [x for x in all_stroke_files if '12months' not in x]
        stroke_files.sort()
        
        #Convert to DataFrame, cull subs that aren't in preprocessed data
        stroke_files = pd.DataFrame(stroke_files, columns = ['files'])
        if n == 0:
            stroke_files['initials'] = os.path.split(os.path.split(stroke_files['files'][n])[0])[0][-2:]
        else:
            stroke_files['initials'] = os.path.split(os.path.split(stroke_files['files'][n])[0])[0][-2:]
        cull_list = stroke_files['initials'].isin(df['initials'])
        df_cull = df[cull_list]
        df_cull.rename(columns = {'Unnamed: 0': 'old_idx'}, inplace = True)
        #Specially remove the second HK (15-28) who refused to die...
        #df_cull = df_cull[df_cull['old_idx'] != 194]
        #df_cull.reset_index(inplace = True, drop = True)
        try:
            np.all(df_cull['initials'] == stroke_files['initials']) == 1       
            stroke_files['madrs_3m'] = df_cull['MADRS_score_3mth']
        except:
            print 'Data are not of identical length. Error in script.'
    
    
    #ADD LOW / HIGH GROUPS TO CULLED DF
    
        stroke_files['group'] = np.zeros(len(stroke_files))
    
        stroke_files['group'][stroke_files['madrs_3m']<=8] = 1
        stroke_files['group'][stroke_files['madrs_3m']>8] = 2
        stroke_files.to_csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_files_madrs_grouped.csv')
    
    #SPLIT GROUPS & MAKE COVARIATE DEFINITIONS
        stroke_low = stroke_files[stroke_files['madrs_3m']<=8]
        stroke_high = stroke_files[stroke_files['madrs_3m']>8]
    
    #CALC MADRS STANDARDIZATION (0 - 1)
        lmin = stroke_low.madrs_3m.min()
        lmax = stroke_low.madrs_3m.max()
        stroke_low['madrs_3m_stand'] = (stroke_low.madrs_3m - lmin) / (lmax - lmin)
    
        hmin = stroke_high.madrs_3m.min()
        hmax = stroke_high.madrs_3m.max()
        stroke_high['madrs_3m_stand'] = (stroke_high.madrs_3m - hmin) / (hmax - hmin)
        
        if stand == 1:
            outdir = outdir[:-1] + '_stand'
            cov1 = np.zeros(len(stroke_files))
            cov1[:len(stroke_low)] = stroke_low['madrs_3m_stand']
            cov1 = {'vector': cov1.tolist(), 'name': 'low_madrs', 'interaction': 1,  'centering': 5} 
            
            cov2 = np.zeros(len(stroke_files))
            cov2[len(stroke_low):] = stroke_high['madrs_3m_stand']
            cov2 = {'vector': cov2.tolist(), 'name': 'high_madrs', 'interaction': 1,  'centering': 5}
        
        elif stand == 0:
            cov1 = np.zeros(len(stroke_files)).astype(int)
            cov1[:len(stroke_low)] = stroke_low['madrs_3m'].astype(int)
            cov1 = {'vector': cov1.tolist(), 'name': 'low_madrs', 'interaction': 1,  'centering': 5} 
            
            cov2 = np.zeros(len(stroke_files)).astype(int)
            cov2[len(stroke_low):] = stroke_high['madrs_3m'].astype(int)
            cov2 = {'vector': cov2.tolist(), 'name': 'high_madrs', 'interaction': 1,  'centering': 5}
                    
        else:
            print('Invalid value %i. Please enter either 0 for non-standardised or 1 for standardised' %(stand),)
            
    
    
    #INITIALISE T-TEST OBJECT
        ttest = Node(TwoSampleTTestDesign(), name = 'TwoSampleT')
        ttest.inputs.group1_files = stroke_low['files'].tolist()
        ttest.inputs.group2_files = stroke_high['files'].tolist()
        ttest.inputs.covariates = [cov1,cov2]
        #ttest.inputs.explicit_mask_file = mask_file
        
    #Estimate model (aka betas)
        modelEst = Node(EstimateModel(), name = 'EstimateModel')
        modelEst.inputs.estimation_method = {'Classical' : 1}
    
    #Estimate contrasts    
        conEst = Node(EstimateContrast(), name = 'EstimateContrasts')
        
        con_1 = ('all','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [1.0, 1.0, 0.0, 0.0])
        con_2 = ('low','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [1.0, 0.0, 0.0, 0.0])
        con_3 = ('high','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [0.0, 1.0, 0.0, 0.0])
        con_4 = ('low GT high','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [1.0, -1.0, 0.0, 0.0])
        con_5 = ('high GT low','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [-1.0, 1.0, 0.0, 0.0])
        con_6 = ('low vs. high','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [0.0, 0.0, 1.0, -1.0])
        con_7 = ('high vs. low','T', ['Group_{1}', 'Group_{2}', 'low_madrs', 'high_madrs'], [0.0, 0.0, -1.0, 1.0])
    
        contrasts = [con_1, con_2, con_3, con_4, con_5, con_6, con_7]
        
        conEst.inputs.contrasts = contrasts
        conEst.inputs.group_contrast = True
        
    #RUN
        l2analysis=Workflow(name = 'l2analysis')
        l2analysis.base_dir = outdir
            
        l2analysis.connect([(ttest,modelEst,[('spm_mat_file','spm_mat_file')]),
                            (modelEst,conEst,[('spm_mat_file','spm_mat_file'),
                                              ('beta_images','beta_images'),
                                              ('residual_image','residual_image')]),
                                              ])
        
        l2analysis.write_graph(graph2use='colored')
        l2analysis.run('MultiProc', plugin_args={'n_procs':1})
    
     


###TWO SAMPLE T###
#Compare controls with prepare stroke (18 vs. 63...)
def prepare_indT_cs():
    outdirs = ['/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_con_vs_stroke/alff_analysis/', '/home/peter/Desktop/prepare/rest/output/sec_level/two_sample_con_vs_stroke/falff_analysis/']
    
    for n,outdir in enumerate(outdirs):
    
        stroke_hdir = '/home/peter/Desktop/prepare/rest/output/'
        all_stroke_files = glob.glob(stroke_hdir + '/*/f_alff/' + maps[n])
        stroke_files = [x for x in all_stroke_files if '12months' not in x]
        
        healthy_hdir = '/home/peter/Desktop/Connect/rest/output/'
        all_healthy_files = glob.glob(healthy_hdir + '/D_H*/f_alff/' + maps[n])
        healthy_files = [x for x in all_healthy_files if 'P2' not in x]
        
        ttest = Node(TwoSampleTTestDesign(), name = 'TwoSampleT')
        ttest.inputs.group1_files = healthy_files
        ttest.inputs.group2_files = stroke_files
        
        modelEst = Node(EstimateModel(), name = 'EstimateModel')
        modelEst.inputs.estimation_method = {'Classical' : 1}
        
        conEst = Node(EstimateContrast(), name = 'EstimateContrasts')
        
        con_1 = ('Controls','T', ['Group_{1}', 'Group_{2}'], [1.0, 0.0])
        con_2 = ('Patients','T', ['Group_{1}','Group_{2}'], [0.0, 1.0])
        con_3 = ('Controls>Patients','T',['Group_{1}','Group_{2}'],[1.0, -1.0])
        con_4 = ('Patients>Controls','T',['Group_{1}','Group_{2}'],[-1.0, 1.0])
        
        contrasts = [con_1, con_2, con_3, con_4]
        
        conEst.inputs.contrasts = contrasts
        conEst.inputs.group_contrast = True
        
        l2analysis=Workflow(name = 'l2analysis')
        l2analysis.base_dir = outdir
            
        l2analysis.connect([(ttest,modelEst,[('spm_mat_file','spm_mat_file')]),
                            (modelEst,conEst,[('spm_mat_file','spm_mat_file'),
                                              ('beta_images','beta_images'),
                                              ('residual_image','residual_image')]),
                                              ])
        
        l2analysis.write_graph(graph2use='colored')
        l2analysis.run('MultiProc', plugin_args={'n_procs':1})
