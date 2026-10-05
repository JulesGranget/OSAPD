


import scipy
import numpy as np





################################
######## GENERAL PARAMS ######## 
################################

srate = 500

sujet_pilot = ['test']

sujet_list = ['P01', 'P02', 'P03', 'P04']


conditions_oc = ['oc300', 'oc600']
conditions_cycle = ['RB', 'O']
conditions_dbs = ['on', 'off']

aux_chanlist = ['resp', 'co2', 'gsm']

freq_band_dict = {'theta' : [4,8], 'alpha' : [8,12], 'beta' : [12,50], 'gamma' : [60, 150]}

ref_type_list = ['mn', 'bi']



########################################
######## PATH DEFINITION ########
########################################

import socket
import os
import platform
 
PC_OS = platform.system()
PC_ID = socket.gethostname()
init_workdir = os.getcwd()

if PC_ID == 'jules-ubuntu1':

    path_main_workdir = '/home/jules/Documents/OSAPD_JULES/Scripts_OSAPD'
    path_general = '/home/jules/Documents/OSAPD_JULES'
    n_core = 25

    
path_data = os.path.join(path_general, 'Data')
path_prep = os.path.join(path_general, 'Analyses', 'preprocessing')
path_anatomy = os.path.join(path_general, 'Data', 'anatomy')
path_precompute = os.path.join(path_general, 'Analyses', 'precompute') 
path_results = os.path.join(path_general, 'Analyses', 'results') 
path_memmap = os.path.join(path_general, 'memmap') 

os.chdir(init_workdir)



################################
######## PREP PARAMS ########
################################ 



prep_step_wb = {
'reref' : {'execute': False, 'params' : ['TP9', 'TP10']}, #chan = chan to reref
'mean_centered' : {'execute': True},
'line_noise_removing' : {'execute': True},
'high_pass' : {'execute': True, 'params' : {'l_freq' : 0.1, 'h_freq': None}},
'low_pass' : {'execute': True, 'params' : {'l_freq' : None, 'h_freq': 200}},
'csd_computation' : {'execute': False},
'ICA_computation' : {'execute': False},
'average_reref' : {'execute': False},
}

bipolarize_chan_i_lim = 8 #non zero based


################################
######## PRECOMPUTE TF ########
################################

#### chunk data
chunk_time = 2 #sec
exclude_cycle_duration_jules = 2 #sec
exclude_inspi_duration_jules = 1 #sec

#### stretch
stretch_point_TF = 250
stretch_TF_auto = False
ratio_stretch_TF = 0.50

#### TF & ITPC
nfrex = 150
ncycle_list = [7, 41]
freq_list = [2, 150]
srate_dw = 10
wavetime = np.arange(-3,3,1/srate)
frex = np.logspace(np.log10(freq_list[0]), np.log10(freq_list[1]), nfrex) 
cycles = np.logspace(np.log10(ncycle_list[0]), np.log10(ncycle_list[1]), nfrex).astype('int')
Pxx_wavelet_norm = 1000

#### STATS
n_surrogates_tf = 1000
tf_percentile_sel_stats = 2 # for both side
tf_stats_percentile_cluster = 95
norm_method = 'rscore'# 'zscore', 'dB'
tf_stats_percentile_cluster_allplot = 99
df_extraction_Cxy = 0.02 #Hz around respi median


#### plot
tf_plot_percentile_scale = 1 #for one side
percentile_plot = [1,99]





################################
######## ERP EXTRACTION ########
################################

#### chunk data
chunk_time_prepost = 5 #sec

#### HFA
HFA_freq = [70, 150]

#### stats
ERP_n_surrogate = 500

