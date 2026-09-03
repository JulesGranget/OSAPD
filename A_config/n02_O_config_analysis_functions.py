

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
import tkinter
import scipy.signal
import mne
import pandas as pd
import sys
import stat
import subprocess
import scipy.stats
import xarray as xr
import physio
import paramiko
import getpass
import cv2
import statsmodels
import seaborn as sns
#install netcdf4

import neurokit2 as nk

from A_config.n01_O_config_params import *
from A_config.n03_O_patient_info import *


debug = False






################################
######## WAVELETS ########
################################


def get_wavelets(srate=500):

    #### compute wavelets
    wavetime = np.arange(-3,3,1/srate)
    wavelets = np.zeros((nfrex, len(wavetime)), dtype=complex)

    # create Morlet wavelet family
    for fi in range(nfrex):
        
        s = cycles[fi] / (2*np.pi*frex[fi])
        gw = np.exp(-wavetime**2/ (2*s**2)) 
        sw = np.exp(1j*(2*np.pi*frex[fi]*wavetime))
        mw =  gw * sw

        wavelets[fi,:] = mw

    if debug:

        plt.plot(np.sum(np.abs(wavelets),axis=1))
        plt.show()

        plt.pcolormesh(np.real(wavelets))
        plt.show()

        plt.plot(np.real(wavelets)[0,:])
        plt.show()

    return wavelets



#freq = freq_band_fc[band]
def get_wavelets_fc(freq):

    #### select wavelet parameters
    wavelets_mask = (frex >= freq[0]) & (frex <= freq[-1])
    frex_list = frex[wavelets_mask]
    ncycle_list = cycles[wavelets_mask]
    nfrex_freq = wavelets_mask.sum()

    if freq[0] < 45:
        wavetime = np.arange(-2,2,1/srate)

    if freq[0] > 45:
        wavetime = np.arange(-.5,.5,1/srate)

    #### compute wavelets
    wavelets = np.zeros((nfrex_freq,len(wavetime)) ,dtype=complex)

    # create Morlet wavelet family
    for fi in range(nfrex_freq):
        
        s = ncycle_list[fi] / (2*np.pi*frex_list[fi])
        gw = np.exp(-wavetime**2/ (2*s**2)) 
        sw = np.exp(1j*(2*np.pi*frex_list[fi]*wavetime))
        mw =  gw * sw

        wavelets[fi,:] = mw

    if debug:

        plt.plot(np.sum(np.abs(wavelets),axis=1))
        plt.show()

        plt.pcolormesh(np.real(wavelets))
        plt.show()

        plt.plot(np.real(wavelets)[0,:])
        plt.show()

    return wavelets





############################
######## LOAD DATA ########
############################




def load_data_sujet(sujet, cond):

    path_source = os.getcwd()
    
    os.chdir(path_prep)

    raw = mne.io.read_raw_fif(f'{sujet}_{cond}.fif', preload=True, verbose='critical')

    data = raw.get_data()

    #### go back to path source
    os.chdir(path_source)

    #### free memory
    del raw

    return data




def load_trial_data_sujet(sujet, trial):

    _path_load_file = os.path.join(path_prep, sujet, 'trial_exports')
    _data = xr.open_dataarray(os.path.join(_path_load_file, f'{sujet}_{trial}.nc'))

    return _data




def get_chanlocalist(sujet):

    path_load_loca = os.path.join(path_prep, 'anatomy', sujet, f"{sujet}_locafile.xlsx")
    df_localist_raw = pd.read_excel(path_load_loca)
    df_localist_sel = df_localist_raw.query(f"select == 1")
    localist = df_localist_sel
    chanlist = df_localist_sel['chan'].values

    return chanlist, localist




def get_alltrials_sujet(sujet):

    path_load_data_trial = os.path.join(path_prep, 'export', sujet)

    data_alltrial = {}
    oc_list_alltrial = {}

    for trial_i, trial_cond in enumerate(trial_list_dict[sujet]):

        data_alltrial[trial_i] = xr.load_dataarray(os.path.join(path_load_data_trial, f"{sujet}_trial{trial_i+1}_preproc.nc"))
        oc_list_alltrial[trial_i] = pd.read_excel(os.path.join(path_load_data_trial, f"{sujet}_trial{trial_i+1}_oc_trig.xlsx"))

    return data_alltrial, oc_list_alltrial






def get_df_loca(sujet):

    path_plotloca = os.path.join(path_prep, sujet, 'anatomy', f"{sujet}_plot_loca.xlsx")
    df_loca = pd.read_excel(path_plotloca)

    return df_loca







########################################
######## LOAD RESPFEATURES ########
########################################


def get_respfeatures(sujet):
    
    path_load_respfeatures = os.path.join(path_precompute, 'respi', 'respfeatures')
    respfeature_allcond = pd.read_excel(os.path.join(path_load_respfeatures, f"{sujet}_respfeatures.xlsx")).drop(columns=['Unnamed: 0'])

    return respfeature_allcond








################################
######## STRETCH ########
################################
    


#resp_features, nb_point_by_cycle = respfeatures[cond], stretch_point_ERP
def stretch_data(resp_features, nb_point_by_cycle, data, srate):

    #### params
    cycle_times = resp_features[['inspi_time', 'expi_time', 'next_inspi_time']].values
    mean_cycle_duration = np.mean(resp_features[['inspi_duration', 'expi_duration']].values, axis=0)
    mean_inspi_ratio = mean_cycle_duration[0]/mean_cycle_duration.sum()
    times = np.arange(0,data.shape[0])/srate

    #### stretch
    if stretch_TF_auto:

        cycles = physio.deform_traces_to_cycle_template(data.reshape(-1,1), times, cycle_times, points_per_cycle=nb_point_by_cycle, 
                segment_ratios=mean_inspi_ratio, output_mode='stacked')
    else:
        
        cycles = physio.deform_traces_to_cycle_template(data.reshape(-1,1), times, cycle_times, points_per_cycle=nb_point_by_cycle, 
                segment_ratios=ratio_stretch_TF, output_mode='stacked')

    #### reshape
    if np.iscomplex(data[0]):
        data_stretch = np.zeros(( cycles.shape[0], nb_point_by_cycle ), dtype='complex')
    else:
        data_stretch = np.zeros(( cycles.shape[0], nb_point_by_cycle ))

    for cycle_i in range(cycles.shape[0]):

        data_stretch[cycle_i, :] = cycles[cycle_i,:].reshape(-1)

    #### inspect
    if debug == True:

        plt.plot(data_stretch.mean(axis=0))
        plt.show()

    return data_stretch, mean_inspi_ratio

    




#resp_features, nb_point_by_cycle, data, srate = respfeature_trial, stretch_point_TF, tf_allconv_norm[chan_i], srate
def stretch_data_tf(resp_features, nb_point_by_cycle, data, srate):

    #### params
    cycle_times = resp_features[['inspi_time', 'expi_time', 'next_inspi_time']].values
    mean_cycle_duration = np.mean(resp_features[['inspi_duration', 'expi_duration']].values, axis=0)
    mean_inspi_ratio = mean_cycle_duration[0]/mean_cycle_duration.sum()
    times = np.arange(0,data.shape[1])/srate

    #### stretch
    if stretch_TF_auto:

        cycles = physio.deform_traces_to_cycle_template(data.T, times, cycle_times, points_per_cycle=nb_point_by_cycle, 
                segment_ratios=mean_inspi_ratio, output_mode='stacked')
    else:
        
        cycles = physio.deform_traces_to_cycle_template(data.T, times, cycle_times, points_per_cycle=nb_point_by_cycle, 
                segment_ratios=ratio_stretch_TF, output_mode='stacked')

    #### reshape
    if np.iscomplex(data[0,0]):
        data_stretch = np.zeros(( cycles.shape[0], data.shape[0], nb_point_by_cycle ), dtype='complex')
    else:
        data_stretch = np.zeros(( cycles.shape[0], data.shape[0], nb_point_by_cycle ))

    for cycle_i in range(cycles.shape[0]):

        data_stretch[cycle_i, :, :] = cycles[cycle_i,:,:].T

    #### inspect
    if debug == True:

        plt.pcolormesh(np.mean(data_stretch, axis=0))
        plt.show()

    return data_stretch, mean_inspi_ratio









########################################
######## CHANGE NAME CSV TRC ########
########################################


def modify_loca_name(loca_list):

    # loca_list = np.concatenate([get_chanlist(sujet)[1] for sujet in sujet_list])

    # np.unique([_loca[:3] for _loca in loca_list])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'ctx'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == '3rd'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'Bra'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'Lef'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'Rig'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'par'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'inf'])

    # np.unique([_loca[7:] for _loca in loca_list if _loca[:3] == 'Ctx'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'ctx'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == '3rd'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'Bra'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'Lef'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'Rig'])
    # np.unique([_loca for _loca in loca_list if _loca[:3] == 'par'])

    loca_list_corrected = []

    for _loca in loca_list:

        if _loca[:3] == '3rd':
            _type = 'csf'
            _loca_corr = _loca
            _side = 'csf'

        elif _loca[:3] == 'Bra':
            _type = 'gm'
            _loca_corr = _loca
            _side = _loca

        elif _loca[:3] == 'ctx':
            _type = 'gm'
            if _loca[4:6] == 'lh':
                _side = 'left'
            elif _loca[4:6] == 'rh':
                _side = 'right'
            
            if _loca[7:].split('-')[-1] in ['ant', 'pos']:
                _loca_corr = _loca[7:]
            else:
                _loca_corr = _loca[7:].split('-')[0]

        elif _loca[:3] == 'Lef':
            _side = 'left'
            if _loca[5:].find('Cerebral-White-Matter') != -1:
                _type = 'wm'
                _loca_corr = 'WM'
            elif _loca[5:].find('Inf-Lat-Vent') != -1:
                _type = 'csf'
                _loca_corr = 'Inf-Lat-Vent'
            else:
                if _loca[5:].split('-')[-1] in ['ant', 'pos']:
                    _type = 'gm'
                    _loca_corr = _loca[5:]
                else:
                    _type = 'gm'
                    _loca_corr = _loca[5:].split('-')[0]

        elif _loca[:3] == 'Rig':
            _side = 'right'
            if _loca[6:].find('Cerebral-White-Matter') != -1:
                _type = 'wm'
                _loca_corr = 'WM'
            elif _loca[6:].find('Inf-Lat-Vent') != -1:
                _type = 'csf'
                _loca_corr = 'Inf-Lat-Vent'
            else:
                if _loca[6:].split('-')[-1] in ['ant', 'pos']:
                    _type = 'gm'
                    _loca_corr = _loca[6:]
                else:
                    _type = 'gm'
                    _loca_corr = _loca[6:].split('-')[0]
            
        elif _loca[:3] == 'Unk':
            _type = 'unknown'
            _loca_corr = 'unknown'
            _side = 'unknown'

        else:
            _type = f'UNSORTED_{_loca}'
            _loca_corr = f'UNSORTED_{_loca}'
            _side = f'UNSORTED_{_loca}'

        loca_list_corrected.append([_loca_corr, _type, _side])

    loca_list_corrected = np.array(loca_list_corrected)

    return loca_list_corrected











########################################
######## MI ANALYSIS FUNCTIONS ########
########################################


#x = data_mat[sujet_i,pair_i,:]
def shuffle_sig(x):

    cut = np.random.randint(low=0, high=len(x), size=1)[0]
    x_cut1 = x[:cut]
    x_cut2 = x[cut:]
    x_shift = np.concatenate((x_cut2, x_cut1), axis=0)

    return x_shift
    

def shuffle_CycleFreq(x):

    cut = int(np.random.randint(low=0, high=len(x), size=1))
    x_cut1 = x[:cut]
    x_cut2 = x[cut:]*-1
    x_shift = np.concatenate((x_cut2, x_cut1), axis=0)

    return x_shift
    

def shuffle_Cxy(x):
   half_size = x.shape[0]//2
   ind = np.random.randint(low=0, high=half_size)
   x_shift = x.copy()
   
   x_shift[ind:ind+half_size] *= -1
   if np.random.rand() >=0.5:
       x_shift *= -1

   return x_shift


def Kullback_Leibler_Distance(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.sum(np.where(a != 0, a * np.log(a / b), 0))

def Shannon_Entropy(a):
    a = np.asarray(a, dtype=float)
    return - np.sum(np.where(a != 0, a * np.log(a), 0))

def Modulation_Index(distrib, show=False, verbose=False):
    distrib = np.asarray(distrib, dtype = float)
    
    if verbose:
        if np.sum(distrib) != 1:
            print(f'(!)  The sum of all bins is not 1 (sum = {round(np.sum(distrib), 2)})  (!)')
        
    N = distrib.size
    uniform_distrib = np.ones(N) * (1/N)
    mi = Kullback_Leibler_Distance(distrib, uniform_distrib) / np.log(N)
    
    if show:
        bin_width_deg = 360 / N
        
        doubled_distrib = np.concatenate([distrib,distrib] )
        x = np.arange(0, doubled_distrib.size*bin_width_deg, bin_width_deg)
        fig, ax = plt.subplots(figsize = (8,4))
        
        doubled_uniform_distrib = np.concatenate([uniform_distrib,uniform_distrib] )
        ax.scatter(x, doubled_uniform_distrib, s=2, color='r')
        
        ax.bar(x=x, height=doubled_distrib, width = bin_width_deg/1.1, align = 'edge')
        ax.set_title(f'Modulation Index = {round(mi, 4)}')
        ax.set_xlabel(f'Phase (Deg)')
        ax.set_ylabel(f'Amplitude (Normalized)')
        ax.set_xticks([0,360,720])

    return mi

def Shannon_MI(a):
    a = np.asarray(a, dtype = float)
    N = a.size
    kl_divergence_shannon = np.log(N) - Shannon_Entropy(a)
    return kl_divergence_shannon / np.log(N)



def get_MVL(x):
    _phase = np.arange(0, x.shape[0])*2*np.pi/x.shape[0]
    complex_vec = x*np.exp(1j*_phase) # ici sous la forme du module * angle, r * phi

    MVL = np.abs(np.mean(complex_vec))
    
    if debug:
        fig = plt.figure()
        ax = fig.add_subplot(projection='polar')
        ax.scatter(complex_vec.real, complex_vec.imag)
        ax.scatter(np.mean(complex_vec.real), np.mean(complex_vec.imag), linewidth=3, color='r')
        plt.show()

    return MVL




def get_MI_2sig(x, y):

    #### Freedman and Diaconis rule
    nbins_x = int(np.ceil((x.max() - x.min()) / (2 * scipy.stats.iqr(x)*(x.size**(-1/3)))))
    nbins_y = int(np.ceil((y.max() - y.min()) / (2 * scipy.stats.iqr(y)*(y.size**(-1/3)))))

    #### compute proba
    hist_x = np.histogram(x,bins = nbins_x)[0]
    hist_x = hist_x/np.sum(hist_x)
    hist_y = np.histogram(y,bins = nbins_y)[0]
    hist_y = hist_y/np.sum(hist_y)

    hist_2d = np.histogram2d(x, y, bins=[nbins_x, nbins_y])[0]
    hist_2d = hist_2d / np.sum(hist_2d)

    #### compute MI
    E_x = 0
    E_y = 0
    E_x_y = 0

    for p in hist_x:
        if p!=0 :
            E_x += -p*np.log2(p)

    for p in hist_y:
        if p!=0 :
            E_y += -p*np.log2(p)

    for p0 in hist_2d:
        for p in p0 :
            if p!=0 :
                E_x_y += -p*np.log2(p)

    MI = E_x+E_y-E_x_y

    return MI



def get_ISPC_2sig(x, y):

    ##### collect "eulerized" phase angle differences
    phase_angle_diff = np.exp(1j*(np.angle(x)-np.angle(y)))

    ##### compute ISPC
    ISPC = np.abs( np.mean(phase_angle_diff) )

    return ISPC



def get_WPLI_2sig(x, y):

    sxy = x * np.conj(y)

    # Extract imaginary part (which is sin(phase difference))
    im_part = np.imag(sxy)
    
    # Compute the weighted phase lag index (wPLI)
    numerator = np.abs(np.mean(im_part))  # Mean of the sign of the imaginary part
    denominator = np.mean(np.abs(im_part))  # Mean of the absolute imaginary part
    
    WPLI = numerator / denominator

    return WPLI


def get_DWPLI_2sig_DEBUG(x, y):

    im_csd = np.imag(x * np.conj(y))  # Imaginary part of the cross-spectral density
    N = len(im_csd)  # Number of samples

    sum_im_csd = np.sum(im_csd) / N
    sum_abs_im_csd = np.sum(np.abs(im_csd)) / N
    sum_sq_im_csd = np.sum(im_csd**2) / N

    num = sum_im_csd**2 - sum_sq_im_csd
    denom = sum_abs_im_csd**2 - sum_sq_im_csd

    # Avoid division by zero
    if denom == 0:
        return 0.0

    DWPLI = num / denom

    return DWPLI


def get_DWPLI_2sig(x, y):

    im_csd = np.imag(x * np.conj(y))  # Imaginary part of the cross-spectral density

    num = np.mean(im_csd)**2 - np.mean(im_csd**2)
    denom = np.mean(np.abs(im_csd))**2 - np.mean(im_csd**2)

    # Avoid division by zero
    if denom == 0:
        return 0.0

    DWPLI = num / denom

    return DWPLI

def get_Cxy_2sig(x, y):

    xy = x * np.conj(y)
    xx = x * np.conj(x)
    yy = y * np.conj(y)

    num = np.mean(np.abs(xy))**2
    denom = np.mean(np.abs(xx)) * np.mean(np.abs(yy))

    # Prevent division by zero
    if denom <= 0:
        return 0.0  

    Cxy = num / denom

    return Cxy






########################################
######## SCRIPT ADVANCEMENT ########
########################################


def print_advancement(i, i_final, steps=[25, 50, 75]):

    steps_i = {}
    for step in steps:

        step_i = 0
        while (step_i/i_final*100) < step:
            step_i += 1

        steps_i[step] = step_i

    for step, step_i in steps_i.items():

        if i == step_i:
            print(f'{step}%', flush=True)






################################
######## NORMALIZATION ########
################################


def zscore(x):

    x_zscore = (x - x.mean()) / x.std()

    return x_zscore




def zscore_mat(x):

    _zscore_mat = (x - x.mean(axis=1).reshape(-1,1)) / x.std(axis=1).reshape(-1,1)

    return _zscore_mat




def rscore(x, median=None, mad=None, axis=0):
    """
    Computes the robust z-score for a vector or 2D matrix along a specified axis.
    
    Parameters:
    - x: np.ndarray, input data (1D or 2D).
    - median: float or np.ndarray, precomputed median (optional).
    - mad: float or np.ndarray, precomputed median absolute deviation (optional).
    - axis: int, axis along which to compute the median and MAD (default: 0).
    
    Returns:
    - rzscore_x: np.ndarray, robust z-score of the input along the given axis.
    """
    if median is None:
        median = np.median(x, axis=axis, keepdims=True)  # Compute median along the given axis
    
    if mad is None:
        mad = np.median(np.abs(x - median), axis=axis, keepdims=True)  # Compute MAD along the given axis

    rzscore_x = (x - median) * 0.6745 / mad  # Compute robust z-score

    return rzscore_x

    



def rscore_mat(x):

    mad = np.median(np.abs(x-np.median(x, axis=1).reshape(-1,1)), axis=1) # median_absolute_deviation

    _rscore_mat = (x-np.median(x, axis=1).reshape(-1,1)) * 0.6745 / mad.reshape(-1,1)

    return _rscore_mat






#tf_conv = tf_median_cycle[nchan, :, :]
def norm_tf(sujet, tf_conv, norm_method):

    path_source = os.getcwd()

    chan_list_sel = []
    chan_list_eeg = []

    if norm_method not in ['rscore', 'zscore']:

        #### load baseline
        os.chdir(os.path.join(path_precompute, sujet, 'baselines'))

        baselines = xr.open_dataarray(f'{sujet}_baselines.nc')

    if norm_method == 'dB':

        for n_chan_i, n_chan in enumerate(chan_list_sel):

            tf_conv[n_chan_i,:,:] = 10*np.log10(tf_conv[n_chan_i,:,:] / baselines.loc[n_chan, :, 'median'].values.reshape(-1,1))

    if norm_method == 'zscore_baseline':

        for n_chan_i, n_chan in enumerate(chan_list_sel):

            tf_conv[n_chan_i,:,:] = (tf_conv[n_chan_i,:,:] - baselines.loc[n_chan,:,'mean'].values.reshape(-1,1)) / baselines.loc[n_chan,:,'std'].values.reshape(-1,1)
                
    if norm_method == 'rscore_baseline':

        for n_chan_i, n_chan in enumerate(chan_list_sel):

            tf_conv[n_chan_i,:,:] = (tf_conv[n_chan_i,:,:] - baselines.loc[n_chan,:,'median'].values.reshape(-1,1)) * 0.6745 / baselines.loc[n_chan,:,'mad'].values.reshape(-1,1)

    if norm_method == 'zscore':

        for n_chan_i, n_chan in enumerate(chan_list_sel):

            tf_conv[n_chan_i,:,:] = zscore_mat(tf_conv[n_chan_i,:,:])
                
    if norm_method == 'rscore':

        for n_chan_i, n_chan in enumerate(chan_list_sel):

            tf_conv[n_chan_i,:,:] = rscore_mat(tf_conv[n_chan_i,:,:])


    #### verify baseline
    if debug:

        nchan = 0
        nchan_name = chan_list_sel[nchan]

        fig, axs = plt.subplots(ncols=2)
        axs[0].set_title('mean std')
        axs[0].plot(baselines.loc[nchan_name,:,'mean'], label='mean')
        axs[0].plot(baselines.loc[nchan_name,:,'std'], label='std')
        axs[0].legend()
        axs[0].set_yscale('log')
        axs[1].set_title('median mad')
        axs[1].plot(baselines.loc[nchan_name,:,'median'], label='median')
        axs[1].plot(baselines.loc[nchan_name,:,'mad'], label='mad')
        axs[1].legend()
        axs[1].set_yscale('log')
        plt.show()

        tf_test = tf_conv[nchan,:,:int(tf_conv.shape[-1]/10)].copy()

        fig, axs = plt.subplots(nrows=6)
        fig.set_figheight(10)
        fig.set_figwidth(15)

        percentile_sel = 0

        vmin = np.percentile(tf_test.reshape(-1),percentile_sel)
        vmax = np.percentile(tf_test.reshape(-1),100-percentile_sel)
        im = axs[0].pcolormesh(tf_test, vmin=vmin, vmax=vmax)
        axs[0].set_title('raw')
        fig.colorbar(im, ax=axs[0])

        tf_baseline = 10*np.log10(tf_test / baselines.loc[chan_list_eeg[nchan], :, 'median'].values.reshape(-1,1))
        vmin = np.percentile(tf_baseline.reshape(-1),percentile_sel)
        vmax = np.percentile(tf_baseline.reshape(-1),100-percentile_sel)
        im = axs[1].pcolormesh(tf_baseline, vmin=vmin, vmax=vmax)
        axs[1].set_title('db')
        fig.colorbar(im, ax=axs[1])

        tf_baseline = (tf_test - baselines.loc[chan_list_eeg[nchan],:,'mean'].values.reshape(-1,1)) / baselines.loc[chan_list_eeg[nchan],:,'std'].values.reshape(-1,1)
        vmin = np.percentile(tf_baseline.reshape(-1),percentile_sel)
        vmax = np.percentile(tf_baseline.reshape(-1),100-percentile_sel)
        im = axs[2].pcolormesh(tf_baseline, vmin=vmin, vmax=vmax)
        axs[2].set_title('zscore')
        fig.colorbar(im, ax=axs[2])

        tf_baseline = (tf_test - baselines.loc[chan_list_eeg[nchan],:,'median'].values.reshape(-1,1)) / baselines.loc[chan_list_eeg[nchan],:,'mad'].values.reshape(-1,1)
        vmin = np.percentile(tf_baseline.reshape(-1),percentile_sel)
        vmax = np.percentile(tf_baseline.reshape(-1),100-percentile_sel)
        im = axs[3].pcolormesh(tf_baseline, vmin=vmin, vmax=vmax)
        axs[3].set_title('rscore')
        fig.colorbar(im, ax=axs[3])

        tf_baseline = zscore_mat(tf_test)
        vmin = np.percentile(tf_baseline.reshape(-1),percentile_sel)
        vmax = np.percentile(tf_baseline.reshape(-1),100-percentile_sel)
        im = axs[4].pcolormesh(tf_baseline, vmin=vmin, vmax=vmax)
        axs[4].set_title('zscore_mat')
        fig.colorbar(im, ax=axs[4])

        tf_baseline = rscore_mat(tf_test)
        vmin = np.percentile(tf_baseline.reshape(-1),percentile_sel)
        vmax = np.percentile(tf_baseline.reshape(-1),100-percentile_sel)
        im = axs[5].pcolormesh(tf_baseline, vmin=vmin, vmax=vmax)
        axs[5].set_title('rscore_mat')
        fig.colorbar(im, ax=axs[5])

        plt.show()

    os.chdir(path_source)

    return tf_conv









########################################
######## HRV ANALYSIS HOMEMADE ########
########################################



#### params
def get_params_hrv_homemade(srate_resample_hrv):
    
    nwind_hrv = int( 128*srate_resample_hrv )
    nfft_hrv = nwind_hrv
    noverlap_hrv = np.round(nwind_hrv/90)
    win_hrv = scipy.signal.windows.hann(nwind_hrv)
    f_RRI = (.1, .5)

    return nwind_hrv, nfft_hrv, noverlap_hrv, win_hrv, f_RRI




#### RRI, IFR
#ecg_i, ecg_cR, srate, srate_resample = ecg_i, ecg_cR, srate, srate_resample_hrv
def get_RRI_IFR(ecg_cR, srate_resample) :

    cR_sec = ecg_cR # cR in sec

    # RRI computation
    RRI = np.diff(cR_sec)
    RRI = np.insert(RRI, 0, np.median(RRI))
    IFR = (1/RRI)

    # interpolate
    f = scipy.interpolate.interp1d(cR_sec, RRI, kind='quadratic', fill_value="extrapolate")
    cR_sec_resample = np.arange(cR_sec[0], cR_sec[-1], 1/srate_resample)
    RRI_resample = f(cR_sec_resample)

    #plt.plot(cR_sec, RRI, label='old')
    #plt.plot(cR_sec_resample, RRI_resample, label='new')
    #plt.legend()
    #plt.show()

    return RRI, RRI_resample, IFR



def get_fig_RRI_IFR(ecg_i, ecg_cR, RRI, IFR, srate, srate_resample):

    cR_sec = ecg_cR # cR in sec
    times = np.arange(0,len(ecg_i))/srate # in sec

    f = scipy.interpolate.interp1d(cR_sec, RRI, kind='quadratic', fill_value="extrapolate")
    cR_sec_resample = np.arange(cR_sec[0], cR_sec[-1], 1/srate_resample)
    RRI_resample = f(cR_sec_resample)

    fig, ax = plt.subplots()
    ax = plt.subplot(411)
    plt.plot(times, ecg_i)
    plt.title('ECG')
    plt.ylabel('a.u.')
    plt.xlabel('s')
    plt.vlines(cR_sec, ymin=min(ecg_i), ymax=max(ecg_i), colors='k')
    plt.subplot(412, sharex=ax)
    plt.plot(cR_sec, RRI)
    plt.title('RRI')
    plt.ylabel('s')
    plt.subplot(413, sharex=ax)
    plt.plot(cR_sec_resample, RRI_resample)
    plt.title('RRI_resampled')
    plt.ylabel('Hz')
    plt.subplot(414, sharex=ax)
    plt.plot(cR_sec, IFR)
    plt.title('IFR')
    plt.ylabel('Hz')
    #plt.show()

    # in this plot one RRI point correspond to the difference value between the precedent RR
    # the first point of RRI is the median for plotting consideration

    return fig




def get_fig_PSD_LF_HF(Pxx, hzPxx, VLF, LF, HF):

    # PLOT
    fig = plt.figure()
    plt.plot(hzPxx,Pxx)
    plt.ylim(0, np.max(Pxx[hzPxx>0.01]))
    plt.xlim([0,.6])
    plt.vlines([VLF, LF, HF], ymin=min(Pxx), ymax=max(Pxx), colors='r')
    #plt.show()
    
    return fig

        
def get_fig_poincarre(RRI):

    RRI_1 = RRI[1:]
    RRI_1 = np.append(RRI_1, RRI[-1]) 

    fig = plt.figure()
    plt.scatter(RRI, RRI_1)
    plt.xlabel('RR (ms)')
    plt.ylabel('RR+1 (ms)')
    plt.title('Poincarré ')
    plt.xlim(.600,1.)
    plt.ylim(.600,1.)

    return fig
    
#### DeltaHR

#RRI, srate_resample, f_RRI, condition = result_struct[keys_result[0]][1], srate_resample, f_RRI, cond 
def get_dHR(RRI_resample, srate_resample, f_RRI):
    
    times = np.arange(0,len(RRI_resample))/srate_resample

        # stairs method
    #RRI_stairs = np.array([])
    #len_cR = len(cR) 
    #for RR in range(len(cR)) :
    #    if RR == 0 :
    #        RRI_i = cR[RR+1]/srate - cR[RR]/srate
    #        RRI_stairs = np.append(RRI_stairs, [RRI_i*1e3 for i in range(int(cR[RR+1]))])
    #    elif RR != 0 and RR != len_cR-1 :
    #        RRI_i = cR[RR+1]/srate - cR[RR]/srate
    #        RRI_stairs = np.append(RRI_stairs, [RRI_i*1e3 for i in range(int(cR[RR+1] - cR[RR]))])
    #    elif RR == len_cR-1 :
    #        RRI_stairs = np.append(RRI_stairs, [RRI_i*1e3 for i in range(int(len(ecg) - cR[RR]))])

    def find_extrema():

        #### define function for HRV

        return

    peaks, troughs = find_extrema(RRI_resample, srate_resample, f_RRI)
    peaks_RRI, troughs_RRI = RRI_resample[peaks], RRI_resample[troughs]
    peaks_troughs = np.stack((peaks_RRI, troughs_RRI), axis=1)

    fig_verif = plt.figure()
    plt.plot(times, RRI_resample)
    plt.vlines(peaks/srate_resample, ymin=min(RRI_resample), ymax=max(RRI_resample), colors='b')
    plt.vlines(troughs/srate_resample, ymin=min(RRI_resample), ymax=max(RRI_resample), colors='r')
    #plt.show()

    dHR = np.diff(peaks_troughs/srate_resample, axis=1)*1e3

    fig_dHR = plt.figure()
    ax = plt.subplot(211)
    plt.plot(times, RRI_resample*1e3)
    plt.title('RRI')
    plt.ylabel('ms')
    plt.subplot(212, sharex=ax)
    plt.plot(troughs/srate_resample, dHR)
    plt.hlines(np.median(dHR), xmin=min(times), xmax=max(times), colors='m', label='median = {:.3f}'.format(np.median(dHR)))
    plt.legend()
    plt.title('dHR')
    plt.ylabel('ms')
    plt.vlines(peaks/srate_resample, ymin=0, ymax=0.01, colors='b')
    plt.vlines(troughs/srate_resample, ymin=0, ymax=0.01, colors='r')
    plt.tight_layout()
    #plt.show()

    return fig_verif, fig_dHR

#ecg_allcond[cond][odor_i], ecg_cR_allcond[cond][odor_i], prms_hrv
def ecg_analysis_homemade(ecg_i, srate, srate_resample_hrv, fig_token=False):

    #### load params
    nwind_hrv, nfft_hrv, noverlap_hrv, win_hrv, f_RRI = get_params_hrv_homemade(srate_resample_hrv)

    #### load cR
    ecg_cR = scipy.signal.find_peaks(ecg_i, distance=srate*0.5)[0]
    ecg_cR = ecg_cR/srate

    #### verif
    if debug:
        times = np.arange(ecg_i.shape[0])/srate
        plt.plot(times, ecg_i)
        plt.vlines(ecg_cR, ymin=np.min(ecg_i) ,ymax=np.max(ecg_i), colors='r')
        plt.show()


    #### initiate metrics names
    res_list = ['HRV_MeanNN', 'HRV_SDNN', 'HRV_RMSSD', 'HRV_pNN50', 'HRV_LF', 'HRV_HF', 'HRV_LFHF', 'HRV_SD1', 'HRV_SD2', 'HRV_S', 'HRV_rCOV', 'HRV_MAD', 'HRV_MEDIAN']

    #### RRI
    RRI, RRI_resample, IFR = get_RRI_IFR(ecg_i, ecg_cR, srate, srate_resample_hrv)

    HRV_MeanNN = np.mean(RRI)
    
    #### PSD
    VLF, LF, HF = .04, .15, .4
    AUC_LF, AUC_HF, LF_HF_ratio, hzPxx, Pxx = get_PSD_LF_HF(RRI_resample, srate_resample_hrv, nwind_hrv, nfft_hrv, noverlap_hrv, win_hrv, VLF, LF, HF)

    #### descriptors
    MeanNN, SDNN, RMSSD, NN50, pNN50, COV, mad, median = get_stats_descriptors(RRI)

    #### poincarré
    SD1, SD2, Tot_HRV = get_poincarre(RRI)

    #### df
    res_tmp = [HRV_MeanNN*1e3, SDNN*1e3, RMSSD, pNN50*100, AUC_LF/10, AUC_HF/10, LF_HF_ratio, SD1*1e3, SD2*1e3, Tot_HRV*1e6, COV, mad*1e3, median*1e3]
    data_df = {}
    for i, dv in enumerate(res_list):
        data_df[dv] = [res_tmp[i]]

    hrv_metrics_homemade = pd.DataFrame(data=data_df)

    #### for figures

    #### dHR
    if fig_token:
        fig_verif, fig_dHR = get_dHR(RRI_resample, srate_resample_hrv, f_RRI)

    #### fig
    if fig_token:
        fig_RRI = get_fig_RRI_IFR(ecg_i, ecg_cR, RRI, IFR, srate, srate_resample_hrv)
        fig_PSD = get_fig_PSD_LF_HF(Pxx, hzPxx, VLF, LF, HF) 
        fig_poincarre = get_fig_poincarre(RRI)

        fig_list = [fig_RRI, fig_PSD, fig_poincarre, fig_verif, fig_dHR]

        plt.close('all')

        return hrv_metrics_homemade, fig_list

    else:

        return hrv_metrics_homemade



def get_hrv_metrics_win(RRI):

    #### initiate metrics names
    res_list = ['HRV_MeanNN', 'HRV_SDNN', 'HRV_RMSSD', 'HRV_pNN50', 'HRV_SD1', 'HRV_SD2', 'HRV_S', 'HRV_COV', 'HRV_MAD', 'HRV_MEDIAN']

    HRV_MeanNN = np.mean(RRI)
    
    #### descriptors
    MeanNN, SDNN, RMSSD, NN50, pNN50, COV, mad, median = get_stats_descriptors(RRI)

    #### poincarré
    SD1, SD2, Tot_HRV = get_poincarre(RRI)

    #### df
    res_tmp = [HRV_MeanNN*1e3, SDNN*1e3, RMSSD, pNN50*100, SD1*1e3, SD2*1e3, Tot_HRV*1e6, COV, mad*1e3, median*1e3]
    data_df = {}
    for i, dv in enumerate(res_list):
        data_df[dv] = [res_tmp[i]]

    hrv_metrics_homemade = pd.DataFrame(data=data_df)

    return hrv_metrics_homemade




def get_PSD_LF_HF(RRI_resample, prms_hrv, VLF, LF, HF):

    srate_resample, nwind, nfft, noverlap, win = prms_hrv['srate_resample_hrv'], prms_hrv['nwind_hrv'], prms_hrv['nfft_hrv'], prms_hrv['noverlap_hrv'], prms_hrv['win_hrv']

    # DETREND
    RRI_detrend = RRI_resample-np.median(RRI_resample)

    # FFT WELCH
    hzPxx, Pxx = scipy.signal.welch(RRI_detrend, fs=srate_resample, window=win, nperseg=nwind, noverlap=noverlap, nfft=nfft)

    AUC_LF = np.trapz(Pxx[(hzPxx>VLF) & (hzPxx<LF)])
    AUC_HF = np.trapz(Pxx[(hzPxx>LF) & (hzPxx<HF)])
    LF_HF_ratio = AUC_LF/AUC_HF

    return AUC_LF, AUC_HF, LF_HF_ratio, hzPxx, Pxx



def get_stats_descriptors(RRI) :

    MeanNN = np.mean(RRI)

    SDNN = np.std(RRI)

    RMSSD = np.sqrt(np.mean((np.diff(RRI)*1e3)**2))

    NN50 = []
    for RR in range(len(RRI)) :
        if RR == len(RRI)-1 :
            continue
        else :
            NN = abs(RRI[RR+1] - RRI[RR])
            NN50.append(NN)

    NN50 = np.array(NN50)*1e3
    pNN50 = np.sum(NN50>50)/len(NN50)

    mad = np.median( np.abs(RRI-np.median(RRI)) )
    COV = mad / np.median(RRI)

    median = np.median(RRI)

    return MeanNN, SDNN, RMSSD, NN50, pNN50, COV, mad, median


def get_poincarre(RRI):
    RRI_1 = RRI[1:]
    RRI_1 = np.append(RRI_1, RRI[-1]) 

    SD1_val = []
    SD2_val = []
    for RR in range(len(RRI)) :
        if RR == len(RRI)-1 :
            continue
        else :
            SD1_val_tmp = (RRI[RR+1] - RRI[RR])/np.sqrt(2)
            SD2_val_tmp = (RRI[RR+1] + RRI[RR])/np.sqrt(2)
            SD1_val.append(SD1_val_tmp)
            SD2_val.append(SD2_val_tmp)

    SD1 = np.std(SD1_val)
    SD2 = np.std(SD2_val)
    Tot_HRV = SD1*SD2*np.pi

    return SD1, SD2, Tot_HRV



def get_hrv_metrics_homemade(cR_time, prms_hrv, analysis_time='5min'):

    #### get RRI
    cR_sec = cR_time/prms_hrv['srate'] # cR in sec

    if analysis_time == '3min':

        cR_sec_mask = (cR_sec >= 60) & (cR_sec <= 240)
        cR_sec = cR_sec[cR_sec_mask] - 60

    RRI = np.diff(cR_sec)
    RRI = np.insert(RRI, 0, np.median(RRI))

    if debug:
        plt.plot(cR_sec, RRI)
        plt.show()
    
    #### descriptors
    MeanNN, SDNN, RMSSD, NN50, pNN50, COV, mad, median = get_stats_descriptors(RRI)

    #### poincarré
    SD1, SD2, Tot_HRV = get_poincarre(RRI)

    #### PSD
    f = scipy.interpolate.interp1d(cR_sec, RRI, kind='quadratic', fill_value="extrapolate")
    cR_sec_resample = np.arange(cR_sec[0], cR_sec[-1], 1/prms_hrv['srate_resample_hrv'])
    RRI_resample = f(cR_sec_resample)

    if debug:
        plt.plot(cR_sec, RRI, label='raw')
        plt.plot(cR_sec_resample, RRI_resample, label='resampled')
        plt.legend()
        plt.show()

    VLF, LF, HF = .04, .15, .4
    AUC_LF, AUC_HF, LF_HF_ratio, hzPxx, Pxx = get_PSD_LF_HF(RRI_resample, prms_hrv, VLF, LF, HF)

    #### df
    res_tmp = {'HRV_MeanNN' : MeanNN*1e3, 'HRV_SDNN' : SDNN*1e3, 'HRV_RMSSD' : RMSSD, 'HRV_pNN50' : pNN50*100, 'HRV_LF' : AUC_LF/10, 'HRV_HF' : AUC_HF/10, 
               'HRV_LFHF' : LF_HF_ratio, 'HRV_SD1' : SD1*1e3, 'HRV_SD2' : SD2*1e3, 'HRV_S' : Tot_HRV*1e6, 'HRV_COV' : COV, 'HRV_MAD' : mad, 'HRV_MEDIAN' : median}
    
    data_df = {}
    for i, dv in enumerate(prms_hrv['metric_list']):
        data_df[dv] = res_tmp[dv]

    hrv_metrics_homemade = pd.DataFrame([data_df])

    return hrv_metrics_homemade









########################################
######## PERMUTATION STATS ######## 
########################################


# data_baseline, data_cond, n_surr = data_Cxy_baseline[:, chan_i], data_Cxy_cond[:, chan_i], n_surrogates_coh
def get_permutation_2groups(data_baseline, data_cond, n_surr, stat_design='within', mode_grouped='median', mode_generate_surr='percentile', percentile_thresh=[0.5, 99.5]):

    if debug:
        count_baseline, _, _ = plt.hist(data_baseline, bins=50, alpha=0.5, label='baseline', color='b')
        count_cond, _, _ = plt.hist(data_cond, bins=50, alpha=0.5, label='cond', color='r')
        plt.vlines([np.median(data_cond)], ymin=0, ymax=count_cond.max(), color='m', linestyles='--')
        plt.vlines([np.median(data_baseline)], ymin=0, ymax=count_baseline.max(), color='c', linestyles='--')
        plt.legend()
        plt.show()

    n_trials_baselines = data_baseline.shape[0]

    data_shuffle = np.concatenate((data_baseline, data_cond), axis=0)
    n_trial_tot = data_shuffle.shape[0]

    if stat_design == 'within':
        if mode_grouped == 'mean':
            obs_distrib = np.mean(data_baseline - data_cond)
        elif mode_grouped == 'median':
            obs_distrib = np.median(data_cond - data_baseline)
    elif stat_design == 'between':
        if mode_grouped == 'mean':
            obs_distrib = np.mean(data_baseline) - np.mean(data_cond)
        elif mode_grouped == 'median':
            obs_distrib = np.median(data_cond) - np.median(data_baseline)

    surr_distrib = np.zeros((n_surr, 2))

    #surr_i = 0
    for surr_i in range(n_surr):

        #### shuffle
        random_sel = np.random.choice(n_trial_tot, size=n_trial_tot, replace=False)
        data_shuffle_baseline = data_shuffle[random_sel[:n_trials_baselines]]
        data_shuffle_cond = data_shuffle[random_sel[n_trials_baselines:]]

        if mode_grouped == 'mean':
            diff_shuffle = data_shuffle_cond.mean() - data_shuffle_baseline.mean()
        elif mode_grouped == 'median':
            diff_shuffle = np.median(data_shuffle_cond) - np.median(data_shuffle_baseline)

        #### generate distrib
        if mode_generate_surr == 'minmax':
            surr_distrib[surr_i, 0], surr_distrib[surr_i, 1] = diff_shuffle.min(), diff_shuffle.max()
        elif mode_generate_surr == 'percentile':
            surr_distrib[surr_i, 0], surr_distrib[surr_i, 1] = np.percentile(diff_shuffle, percentile_thresh[0]), np.percentile(diff_shuffle, percentile_thresh[1])    

    if debug:
        count, _, _ = plt.hist(surr_distrib[:,0], bins=50, color='k', alpha=0.5)
        count, _, _ = plt.hist(surr_distrib[:,1], bins=50, color='k', alpha=0.5)
        plt.vlines([obs_distrib], ymin=0, ymax=count.max(), label='obs', colors='g')

        plt.vlines([np.percentile(surr_distrib[:,0], 0.5)], ymin=0, ymax=count.max(), label='perc_05_995', colors='r', linestyles='--')
        plt.vlines([np.percentile(surr_distrib[:,1], 99.5)], ymin=0, ymax=count.max(), colors='r', linestyles='--')
        plt.vlines([np.percentile(surr_distrib[:,0], 2.5)], ymin=0, ymax=count.max(), label='perc_025_975', colors='r', linestyles='-.')
        plt.vlines([np.percentile(surr_distrib[:,1], 97.5)], ymin=0, ymax=count.max(), colors='r', linestyles='-.')
        plt.legend()
        plt.show()

    #### thresh
    # surr_dw, surr_up = np.percentile(surr_distrib[:,0], 2.5, axis=0), np.percentile(surr_distrib[:,1], 97.5, axis=0)
    # surr_dw, surr_up = np.percentile(surr_distrib[:,0], 0.5, axis=0), np.percentile(surr_distrib[:,1], 99.5, axis=0)
    surr_dw, surr_up = np.percentile(surr_distrib[:,0], percentile_thresh[0], axis=0), np.percentile(surr_distrib[:,1], percentile_thresh[1], axis=0)

    if obs_distrib < surr_dw or obs_distrib > surr_up:
        stats_res = True
    else:
        stats_res = False

    return stats_res





# data_baseline, data_cond, n_surr = data_baseline_rscore, data_cond_rscore, n_surr_fc
def get_permutation_cluster_1d(data_baseline, data_cond, n_surr, stat_design='within', mode_grouped='median', mode_generate_surr='percentile_time', 
                               mode_select_thresh='percentile_time', percentile_thresh=[2.5, 97.5], size_thresh_alpha=0.01):

    n_trials_baselines = data_baseline.shape[0]
    len_sig = data_baseline.shape[-1]

    data_shuffle = np.concatenate((data_baseline, data_cond), axis=0)
    n_trial_tot = data_shuffle.shape[0]

    if stat_design == 'within':
        if mode_grouped == 'mean':
            obs_distrib = np.mean(data_baseline - data_cond, axis=0)
        elif mode_grouped == 'median':
            obs_distrib = np.median(data_cond - data_baseline, axis=0)
    elif stat_design == 'between':
        if mode_grouped == 'mean':
            obs_distrib = np.mean(data_baseline, axis=0) - np.mean(data_cond, axis=0)
        elif mode_grouped == 'median':
            obs_distrib = np.median(data_cond, axis=0) - np.median(data_baseline, axis=0)

    if mode_generate_surr in ['minmax', 'percentile']:
        surr_distrib = np.zeros((n_surr, 2))
    elif mode_generate_surr == 'percentile_time':
        surr_distrib = np.zeros((n_surr, len_sig))

    if debug:

        if mode_grouped == 'mean':
            data_baseline_grouped = np.mean(data_baseline, axis=0)
            data_cond_grouped = np.mean(data_cond, axis=0)
        elif mode_grouped == 'median':
            data_baseline_grouped = np.median(data_baseline, axis=0)
            data_cond_grouped = np.median(data_cond, axis=0)

        time = np.arange(len_sig)
        rsem_baseline = scipy.stats.median_abs_deviation(data_baseline, axis=0)/np.sqrt(data_baseline.shape[0])
        rsem_cond = scipy.stats.median_abs_deviation(data_cond, axis=0)/np.sqrt(data_cond.shape[0])

        plt.plot(time, data_baseline_grouped, label='baseline', color='c')
        plt.fill_between(time, data_baseline_grouped-rsem_baseline, data_baseline_grouped+rsem_baseline, color='c', alpha=0.5)
        plt.plot(time, data_cond_grouped, label='cond', color='g')
        plt.fill_between(time, data_cond_grouped-rsem_cond, data_cond_grouped+rsem_cond, color='g', alpha=0.5)
        plt.legend()
        plt.show()

    #surr_i = 0
    for surr_i in range(n_surr):

        #### shuffle
        random_sel = np.random.choice(n_trial_tot, size=n_trial_tot, replace=False)
        data_shuffle_baseline = data_shuffle[random_sel[:n_trials_baselines]]
        data_shuffle_cond = data_shuffle[random_sel[n_trials_baselines:]]

        if mode_grouped == 'mean':
            diff_shuffle = np.mean(data_shuffle_cond, axis=0) - np.mean(data_shuffle_baseline, axis=0)
        elif mode_grouped == 'median':
            diff_shuffle = np.median(data_shuffle_cond, axis=0) - np.median(data_shuffle_baseline, axis=0)

        if debug:
            plt.plot(np.mean(data_shuffle_baseline, axis=0), label='baseline')
            plt.plot(np.mean(data_shuffle_cond, axis=0), label='cond')
            plt.legend()
            plt.show()

            plt.hist(np.median(data_shuffle_baseline, axis=0), bins=50, label='baseline', alpha=0.5)
            plt.hist(np.median(data_shuffle_cond, axis=0), bins=50, label='cond', alpha=0.5)
            plt.legend()
            plt.show()

        #### generate distrib
        if mode_generate_surr == 'minmax':
            surr_distrib[surr_i, 0], surr_distrib[surr_i, 1] = diff_shuffle.min(), diff_shuffle.max()
        elif mode_generate_surr == 'percentile':
            surr_distrib[surr_i, 0], surr_distrib[surr_i, 1] = np.percentile(diff_shuffle, 1), np.percentile(diff_shuffle, 99)    
        elif mode_generate_surr == 'percentile_time':
            surr_distrib[surr_i, :] = diff_shuffle

    if debug:
        count, _, _ = plt.hist(surr_distrib[:,0], bins=50, color='k', alpha=0.5)
        count, _, _ = plt.hist(surr_distrib[:,1], bins=50, color='k', alpha=0.5)
        count, _, _ = plt.hist(obs_distrib, bins=50, label='obs', color='g')
        plt.vlines([np.median(surr_distrib[:,0])], ymin=0, ymax=count.max(), label='median', colors='r')
        plt.vlines([np.median(surr_distrib[:,1])], ymin=0, ymax=count.max(), colors='r')
        plt.vlines([np.mean(surr_distrib[:,0])], ymin=0, ymax=count.max(), label='mean', colors='b')
        plt.vlines([np.mean(surr_distrib[:,1])], ymin=0, ymax=count.max(), colors='b')
        plt.vlines([np.percentile(surr_distrib[:,0], 1)], ymin=0, ymax=count.max(), label='perc_1_99', colors='r', linestyles='--')
        plt.vlines([np.percentile(surr_distrib[:,1], 99)], ymin=0, ymax=count.max(), colors='r', linestyles='--')
        plt.vlines([np.percentile(surr_distrib[:,0], 2.5)], ymin=0, ymax=count.max(), label='perc_025_975', colors='r', linestyles='-.')
        plt.vlines([np.percentile(surr_distrib[:,1], 97.5)], ymin=0, ymax=count.max(), colors='r', linestyles='-.')
        plt.legend()
        plt.show()

        plt.plot(obs_distrib)
        plt.hlines([np.median(surr_distrib[:,0])], xmin=0, xmax=len_sig, label='median', colors='r')
        plt.hlines([np.median(surr_distrib[:,1])], xmin=0, xmax=len_sig, colors='r')
        plt.hlines([np.mean(surr_distrib[:,0])], xmin=0, xmax=len_sig, label='mean', colors='b')
        plt.hlines([np.mean(surr_distrib[:,1])], xmin=0, xmax=len_sig, colors='b')
        plt.hlines([np.percentile(surr_distrib[:,0], 0.5)], xmin=0, xmax=len_sig, label='perc_005_995', colors='r', linestyles='--')
        plt.hlines([np.percentile(surr_distrib[:,1], 99.5)], xmin=0, xmax=len_sig, colors='r', linestyles='--')
        plt.hlines([np.percentile(surr_distrib[:,0], 2.5)], xmin=0, xmax=len_sig, label='perc_025_975', colors='r', linestyles='-.')
        plt.hlines([np.percentile(surr_distrib[:,1], 97.5)], xmin=0, xmax=len_sig, colors='r', linestyles='-.')
        plt.hlines([np.percentile(surr_distrib[:,0], 2.5)], xmin=0, xmax=len_sig, label='perc_025_975', colors='r', linestyles='-.')
        plt.hlines([np.percentile(surr_distrib[:,1], 97.5)], xmin=0, xmax=len_sig, colors='r', linestyles='-.')
        plt.legend()
        plt.show()

        plt.plot(obs_distrib)
        plt.plot(np.percentile(surr_distrib, 0.5, axis=0), color='r', linestyle='--')
        plt.plot(np.percentile(surr_distrib, 99.5, axis=0), color='r', linestyle='--')
        plt.plot(np.percentile(surr_distrib, 2.5, axis=0), color='m', linestyle='-.')
        plt.plot(np.percentile(surr_distrib, 97.5, axis=0), color='m', linestyle='-.')
        plt.legend()
        plt.show()

    if mode_select_thresh == 'percentile':
        # surr_dw, surr_up = np.percentile(surr_distrib[:,0], 2.5, axis=0), np.percentile(surr_distrib[:,1], 97.5, axis=0)
        # surr_dw, surr_up = np.percentile(surr_distrib[:,0], 1, axis=0), np.percentile(surr_distrib[:,1], 99, axis=0)
        surr_dw, surr_up = np.percentile(surr_distrib[:,0], percentile_thresh[0], axis=0), np.percentile(surr_distrib[:,1], percentile_thresh[1], axis=0)
    elif mode_select_thresh == 'mean':
        surr_dw, surr_up = np.mean(surr_distrib[:,0], axis=0), np.median(surr_distrib[:,1], axis=0)
    elif mode_select_thresh == 'median':
        surr_dw, surr_up = np.median(surr_distrib[:,0], axis=0), np.median(surr_distrib[:,1], axis=0)
    elif mode_select_thresh == 'percentile_time':
        # surr_dw, surr_up = np.percentile(surr_distrib, 0.5, axis=0), np.percentile(surr_distrib, 99.5, axis=0)
        # surr_dw, surr_up = np.percentile(surr_distrib, 2.5, axis=0), np.percentile(surr_distrib, 97.5, axis=0)
        surr_dw, surr_up = np.percentile(surr_distrib, percentile_thresh[0], axis=0), np.percentile(surr_distrib, percentile_thresh[1], axis=0)

    #### thresh data
    mask = (obs_distrib < surr_dw) | (obs_distrib > surr_up)

    if debug:

        plt.scatter(range(mask.size), mask)
        plt.show()

    if mask.sum() != 0:
    
        #### thresh cluster
        mask_thresh = mask.astype('uint8')
        nb_blobs, im_with_separated_blobs, stats, _ = cv2.connectedComponentsWithStats(mask_thresh)
        #### nb_blobs, im_with_separated_blobs, stats = nb clusters, clusters image with labeled clusters, info on clusters
        sizes = stats[1:, -1]
        nb_blobs -= 1
        # min_size = np.percentile(sizes,size_thresh)  
        min_size = len_sig*size_thresh_alpha  

        if debug:

            count, _, _ = plt.hist(sizes, bins=50, cumulative=True)
            plt.vlines(min_size, ymin=0, ymax=count.max(), colors='r')
            plt.show()

        mask_thresh = np.zeros_like(im_with_separated_blobs)
        for blob in range(nb_blobs):
            if sizes[blob] >= min_size:
                mask_thresh[im_with_separated_blobs == blob + 1] = 1

        mask_thresh = mask_thresh.reshape(-1)

        if debug:

            time = np.arange(data_baseline.shape[-1])
            sem_baseline = data_baseline.std(axis=0)/np.sqrt(data_baseline.shape[0])
            sem_cond = data_cond.std(axis=0)/np.sqrt(data_cond.shape[0])

            plt.plot(time, data_baseline_grouped, label='baseline', color='c')
            plt.fill_between(time, data_baseline_grouped-sem_baseline, data_baseline_grouped+sem_baseline, color='c', alpha=0.5)
            plt.plot(time, data_cond_grouped, label='cond', color='g')
            plt.fill_between(time, data_cond_grouped-sem_cond, data_cond_grouped+sem_cond, color='g', alpha=0.5)
            plt.fill_between(time, data_baseline_grouped.min(), data_cond_grouped.max(), where=mask, color='r', alpha=0.5)
            plt.title('mask not threshed')
            plt.legend()
            plt.show()

            plt.plot(time, data_baseline_grouped, label='baseline', color='c')
            plt.fill_between(time, data_baseline_grouped-sem_baseline, data_baseline_grouped+sem_baseline, color='c', alpha=0.5)
            plt.plot(time, data_cond_grouped, label='cond', color='g')
            plt.fill_between(time, data_cond_grouped-sem_cond, data_cond_grouped+sem_cond, color='g', alpha=0.5)
            plt.fill_between(time, data_baseline_grouped.min(), data_cond_grouped.max(), where=mask_thresh, color='r', alpha=0.5)
            plt.title('mask threshed')
            plt.legend()
            plt.show()

    else:

        mask_thresh = mask

    return mask_thresh




# data_baseline, data_cond, n_surr = tf_stretch_baseline_allsujet, tf_stretch_cond_allsujet, 1000
def get_permutation_cluster_2d(data_baseline, data_cond, n_surr, stat_design='within', mode_grouped='median', mode_generate_surr='percentile_time', 
                               mode_select_thresh='percentile_time', percentile_thresh=[0.5, 99.5], size_thresh_alpha=0.01):

    #### define ncycle
    n_trial_baselines = data_baseline.shape[0]
    n_trial_cond = data_cond.shape[0]
    n_trial_tot = n_trial_baselines + n_trial_cond
    len_sig = data_baseline.shape[-1]

    data_shuffle = np.concatenate((data_baseline, data_cond), axis=0)

    if stat_design == 'within':
        if mode_grouped == 'mean':
            obs_distrib = np.mean(data_cond - data_baseline, axis=0)
        elif mode_grouped == 'median':
            obs_distrib = np.median(data_cond - data_baseline, axis=0)
    elif stat_design == 'between':
        if mode_grouped == 'mean':
            obs_distrib = np.mean(data_cond, axis=0) - np.mean(data_baseline, axis=0)
        elif mode_grouped == 'median':
            obs_distrib = np.median(data_cond, axis=0) - np.median(data_baseline, axis=0)

    if debug:

        plt.pcolormesh(np.median(data_baseline, axis=0))
        plt.show()

        plt.pcolormesh(np.median(data_cond, axis=0))
        plt.show()

        plt.pcolormesh(obs_distrib)
        plt.show()

    #### space allocation
    if mode_generate_surr in ['minmax', 'percentile']:
        surr_distrib = np.zeros((n_surr, 2))
    elif mode_generate_surr == 'percentile_time':
        surr_distrib = np.zeros((n_surr, data_baseline.shape[1], len_sig))

    #surr_i = 0
    for surr_i in range(n_surr):

        print_advancement(surr_i, n_surr, steps=[25, 50, 75])

        #### shuffle
        random_sel = np.random.choice(n_trial_tot, size=n_trial_tot, replace=False)
        data_shuffle_baseline = data_shuffle[random_sel[:n_trial_baselines]]
        data_shuffle_cond = data_shuffle[random_sel[n_trial_baselines:]]

        if stat_design == 'within':
            if mode_grouped == 'mean':
                diff_shuffle = np.mean(data_shuffle_cond - data_shuffle_baseline, axis=0)
            elif mode_grouped == 'median':
                diff_shuffle = np.median(data_shuffle_cond - data_shuffle_baseline, axis=0)
        elif stat_design == 'between':
            if mode_grouped == 'mean':
                diff_shuffle = np.mean(data_shuffle_cond, axis=0) - np.mean(data_shuffle_baseline, axis=0)
            elif mode_grouped == 'median':
                diff_shuffle = np.median(data_shuffle_cond, axis=0) - np.median(data_shuffle_baseline, axis=0)

        if debug:
            plt.pcolormesh(diff_shuffle)
            plt.show()

        #### generate distrib
        if mode_generate_surr == 'minmax':
            surr_distrib[:, surr_i, 0], surr_distrib[:, surr_i, 1] = diff_shuffle.min(axis=1), diff_shuffle.max(axis=1)
        elif mode_generate_surr == 'percentile_time':
            surr_distrib[surr_i] = diff_shuffle

    if mode_select_thresh == 'percentile':
        # surr_dw, surr_up = np.percentile(surr_distrib[:,:,0], 2.5, axis=1), np.percentile(surr_distrib[:,:,1], 97.5, axis=1)
        surr_dw, surr_up = np.percentile(surr_distrib[:,:,0], 1, axis=1), np.percentile(surr_distrib[:,:,1], 99, axis=1)
    elif mode_select_thresh == 'mean':
        surr_dw, surr_up = np.mean(surr_distrib[:,:,0], axis=1), np.median(surr_distrib[:,:,1], axis=1)
    elif mode_select_thresh == 'median':
        surr_dw, surr_up = np.median(surr_distrib[:,:,0], axis=1), np.median(surr_distrib[:,:,1], axis=1)
    elif mode_select_thresh == 'percentile_time':
        surr_dw, surr_up = np.percentile(surr_distrib, percentile_thresh[0], axis=0), np.percentile(surr_distrib, percentile_thresh[1], axis=0)

    if debug:

        bins=50
        counts = np.zeros((obs_distrib.shape[0], bins))
        values = np.zeros((obs_distrib.shape[0], bins+1))
        for row_i in range(obs_distrib.shape[0]):
            counts[row_i,:], values[row_i,:], _ = plt.hist(obs_distrib[row_i,:], bins=bins)
        plt.close('all')

        fig, ax = plt.subplots(figsize=(8, 6))

        X, Y = np.meshgrid(values[0, :-1], np.arange(obs_distrib.shape[0]))  # Mesh grid for pcolormesh

        c = ax.pcolormesh(X, Y, counts, cmap='viridis', shading='auto')

        ax.plot(surr_dw, np.arange(150), color='red', linewidth=2, label="surr_dw")
        ax.plot(surr_up, np.arange(150), color='blue', linewidth=2, label="surr_up")

        ax.set_xlabel("Value Distribution")
        ax.set_ylabel("150 Points")
        ax.set_title("Distribution of Values with Vector Overlays")
        ax.legend()

        fig.colorbar(c, ax=ax, label="Density")

        plt.show()

    #### thresh data
    mask = (obs_distrib < surr_dw) | (obs_distrib > surr_up)

    if debug:

        plt.pcolormesh(mask)
        plt.show()

    if mask.sum() != 0:
    
        #### thresh cluster
        mask_thresh = mask.astype('uint8')
        nb_blobs, im_with_separated_blobs, stats, _ = cv2.connectedComponentsWithStats(mask_thresh)
        #### nb_blobs, im_with_separated_blobs, stats = nb clusters, clusters image with labeled clusters, info on clusters
        sizes = stats[1:, -1]
        nb_blobs -= 1
        # min_size = np.percentile(sizes,size_thresh)  
        min_size = len_sig*size_thresh_alpha  

        if debug:

            count, _, _ = plt.hist(sizes, bins=50, cumulative=True)
            plt.vlines(min_size, ymin=0, ymax=count.max(), colors='r')
            plt.show()

        mask_thresh = np.zeros_like(im_with_separated_blobs)
        for blob in range(nb_blobs):
            blob_vec = np.sum(im_with_separated_blobs == blob + 1, axis=0)
            if np.sum(blob_vec != 0) >= min_size:
                mask_thresh[im_with_separated_blobs == blob + 1] = 1

        if debug:

            fig, ax = plt.subplots()

            ax.pcolormesh(obs_distrib, shading='gouraud', cmap=plt.get_cmap('seismic'))
            ax.contour(mask_thresh, levels=0, colors='g')

            plt.show()

    else:

        mask_thresh = mask

    return mask_thresh






########################
######## FILTER ########
########################

#sig = data
def iirfilt(sig, srate, lowcut=None, highcut=None, order=4, ftype='butter', verbose=False, show=False, axis=0):

    if len(sig.shape) == 1:

        axis = 0

    if lowcut is None and not highcut is None:
        btype = 'lowpass'
        cut = highcut

    if not lowcut is None and highcut is None:
        btype = 'highpass'
        cut = lowcut

    if not lowcut is None and not highcut is None:
        btype = 'bandpass'

    if btype in ('bandpass', 'bandstop'):
        band = [lowcut, highcut]
        assert len(band) == 2
        Wn = [e / srate * 2 for e in band]
    else:
        Wn = float(cut) / srate * 2

    filter_mode = 'sos'
    sos = scipy.signal.iirfilter(order, Wn, analog=False, btype=btype, ftype=ftype, output=filter_mode)

    filtered_sig = scipy.signal.sosfiltfilt(sos, sig, axis=axis)

    return filtered_sig




