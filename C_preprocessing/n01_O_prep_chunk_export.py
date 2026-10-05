
import os
import neo
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
import joblib
import physio
import math
import mne
import scipy

from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *


debug = False








############################################
######## OC EXTRACTION FUNCTION ########
############################################



def get_oc_values(resp_cycle_sig, expi_index, disp_factor=0.375):

    mid_inspi = int(expi_index/2)
    start_inspi_dec, start_inspi_asc = int(mid_inspi-expi_index*disp_factor), int(mid_inspi+expi_index*disp_factor)
    median_inspi_trough = np.median(resp_cycle_sig[start_inspi_dec:start_inspi_asc])
    oc_top_i = np.argmin(resp_cycle_sig[:expi_index])
    oc_top_val = resp_cycle_sig[oc_top_i]

    _oc_ratio = (oc_top_val - median_inspi_trough) / median_inspi_trough
    
    return _oc_ratio, oc_top_val





########################################
######## DATA EXTRACTION ######## 
########################################

#resp, ieeg, resp_trig, ieeg_trig, resp_srate, ieeg_srate, target_srate = resp_chunk, data_ieeg_sync, trigger_SL_trial, data_trig_ieeg_sync, srate_SL, srate_ieeg, srate
def sync_resp_ieeg(resp, ieeg, resp_trig, ieeg_trig, resp_srate, ieeg_srate, target_srate):

    #### clock main vec
    a, b = np.polyfit(resp_trig, ieeg_trig, 1)

    #### transform time_vec
    t_resp = np.arange(resp.shape[-1]) / resp_srate
    t_ieeg = np.arange(ieeg.shape[-1]) / ieeg_srate

    t_resp_sync = a * t_resp + b

    #### syncho strat and stop
    t_start = max(t_resp_sync[0], t_ieeg[0])
    t_stop = min(t_resp_sync[-1], t_ieeg[-1])

    #### common target srate time vector
    t_common = np.arange(t_start, t_stop, 1 / target_srate)

    #### interpolate respiration
    f_resp = scipy.interpolate.interp1d(t_resp_sync, resp, kind="linear", bounds_error=False, fill_value=np.nan)
    resp_resampled = f_resp(t_common)

    if ieeg.shape[0] < 2:

        f_ieeg = scipy.interpolate.interp1d(t_ieeg, ieeg, axis=-1, kind="linear", bounds_error=False, fill_value=np.nan)
        ieeg_resampled = f_ieeg(t_common)

    else:

        ieeg_resampled = []

        for chan_i in range(ieeg.shape[0]):

            f_ieeg = scipy.interpolate.interp1d(t_ieeg, ieeg[chan_i], axis=-1, kind="linear", bounds_error=False, fill_value=np.nan)
            ieeg_resampled.append(f_ieeg(t_common))

        ieeg_resampled = np.array(ieeg_resampled)

    return resp_resampled, ieeg_resampled






#sujet = sujet_list[1]
def extract_data(sujet):

    #### config
    path_data_raw = os.path.join(path_data, sujet, 'neural')
    path_export_resp_fig = os.path.join(path_results, 'respi', 'detection')
    path_resp_data = os.path.join(path_data, sujet, 'physio', f"OSAPD_{sujet}.csv")

    srate = 500
    srate_SL = 100

    #### channel name load
    nchan_extraction = 16
    chanlist_extraction = []

    for chan_i in range(nchan_extraction):
        if chan_i < 9:
            chanlist_extraction.append(f"CECOG_LF_1___0{chan_i+1}")
        else:
            chanlist_extraction.append(f"CECOG_LF_1___{chan_i+1}")

    if sujet in ['P01', 'P02']:
        chan_trig_extraction = ['CDIG_IN_1_Up', 'CDIG_IN_1_Down']
        srate_trig_extraction = 'CDIG_IN_1_KHz'
    else:
        chan_trig_extraction = 'CANALOG_IN_1'
        srate_trig_extraction = 'CANALOG_IN_1_KHz'

    chanlist = chanlist_extraction

    #### breathing load
    df_SL_raw = pd.read_csv(path_resp_data, low_memory=False)
    resp_raw = df_SL_raw['Module 2: Raw Pressure'].values[1:].astype('float')
    trig_SL = df_SL_raw['Digital I/O: Input 2'].values[1:].astype('float')

    if debug:

        plt.plot(scipy.stats.zscore(resp_raw))
        plt.plot(scipy.stats.zscore(trig_SL))
        plt.show()

        plt.plot(df_SL_raw['Serial Port 1: CO2 Concentration'].values[1:].astype('float'))
        plt.plot(df_SL_raw['Digital I/O: Input 1'].values[1:].astype('float'))
        plt.plot(df_SL_raw['Digital I/O: Input 2'].values[1:].astype('float'))
        # plt.plot(scipy.stats.zscore(df_SL_raw['Module 2: Raw Pressure'].values[1:].astype('float')))
        plt.plot(scipy.stats.zscore(df_SL_raw['Module 1: Raw Flow'].values[1:].astype('float')))
        plt.plot(scipy.stats.zscore(df_SL_raw['Module 4: Raw Flow'].values[1:].astype('float')))
        plt.show()

    #### clean trigger
    if debug:

        trigger_spot = np.array([[57900, 58100], [131300, 131350]])
        plt.plot(scipy.stats.zscore(resp_raw))
        plt.plot(scipy.stats.zscore(trig_SL))
        plt.vlines(trigger_spot.reshape(-1), ymin=scipy.stats.zscore(resp_raw).min(), ymax=scipy.stats.zscore(resp_raw).max(), color='r')
        plt.show()

    if sujet in correction_trig_all_patient:
        
        for start, stop in correction_trig_all_patient[sujet]:
            trig_SL[start:stop] = 1

    #### trigger correction
    trig_diff = np.diff(trig_SL)
    trigger_SL = np.where(trig_diff != 0)[0]

    if sujet == 'P01':
        trigger_SL_trial_chunk = np.array([[11822, 42448], [66579, 98191], [108058, 139117], [149954, 179021]])
    else:
        diff_trig_i = np.diff(trigger_SL)
        diff_trig_i = diff_trig_i[diff_trig_i != 1]
        thresh_trial_start_stop = (np.log(diff_trig_i)*-1).max()*2.5
        _trig_i_list_fast = np.where(thresh_trial_start_stop < np.log(np.diff(trigger_SL))*-1)[0]
        diff_fast_trig = np.diff(trigger_SL[_trig_i_list_fast])
        trigger_SL_trial_chunk_i = np.where((diff_fast_trig.max() - diff_fast_trig.min())/10 < diff_fast_trig)[0]
        trigger_SL_trial_chunk_i = np.concat([trigger_SL_trial_chunk_i, [diff_fast_trig.size]])
        trigger_SL_trial_chunk = trigger_SL[_trig_i_list_fast[trigger_SL_trial_chunk_i]+1]
        trigger_SL_trial_chunk = trigger_SL_trial_chunk.reshape(4,-1)

    if debug:
            
        diff_diff_trig = np.log(np.diff(trigger_SL))*-1
        plt.plot(diff_diff_trig)
        plt.hlines(thresh_trial_start_stop, xmin=0, xmax=diff_diff_trig.size, color='r')
        plt.show()

        plt.plot(trig_diff)
        plt.vlines(trigger_SL, ymin=-2, ymax=2, color='r')
        plt.vlines(trigger_SL_trial_chunk.reshape(-1), ymin=-2, ymax=2, color='g')
        plt.show()

    #### extract ieeg
    data_ieeg_alltrial_raw = []
    resp_alltrial = []
    oc_trigger_alltrial = []
    respfeature_alltrial = []

    #trial_i, trial_cond = 1, trial_list_dict[sujet][1]
    for trial_i, trial_cond in enumerate(trial_list_dict[sujet]):

        print(trial_i, trial_cond)
        
        filename = [_file for _file in os.listdir(path_data_raw) if _file[:2] == f"0{trial_i+1}"][0]
        data_raw = scipy.io.loadmat(os.path.join(path_data_raw, filename))
        data_raw = {key.replace('\x00', 'x00'): value for key, value in data_raw.items()}

        if debug:

            data_raw.keys()

            key_open_sig = [_key for _key in data_raw.keys() if _key.find('x00') != -1]
            
            for _key_sig in key_open_sig:
                plt.plot(data_raw[_key_sig].reshape(-1))
                plt.title(_key_sig)
                plt.show()

        data_ieeg = np.array([data_raw[f"{chan}x00"].reshape(-1) for chan in chanlist_extraction])
        srate_ieeg_list = np.array([data_raw[f"{chan}_KHz"].reshape(-1) for chan in chanlist_extraction])

        if ~(srate_ieeg_list == srate_ieeg_list[0]).all():
            raise ValueError('!!! srate different !!!')

        srate_ieeg = int(srate_ieeg_list[0,0]*1e3)

        if debug:

            chan_i = 2

            for chan_i, chan in enumerate(chanlist):

                x = data_ieeg[chan_i]

                frex, Pxx = scipy.signal.welch(x, fs=srate, nperseg=int(srate * 4))
                plt.semilogy(frex, Pxx, label=chan)

            plt.legend()
            plt.show()

        #### ieeg trig
        if sujet in ['P01', 'P02']:

            trial_start = np.array([data_raw[f"{chan}_TimeBegin"].reshape(-1) for chan in chanlist_extraction])
            trial_stop = np.array([data_raw[f"{chan}_TimeEnd"].reshape(-1) for chan in chanlist_extraction])
    
            # (trial_start == trial_start[0]).all() or ~(trial_stop == trial_stop[0]).all()
    
            trial_start = trial_start[0]
            trial_stop = trial_stop[0] - trial_start[0]

            data_trig_ieeg_srate = int(data_raw[f"{srate_trig_extraction}x00"][0,0]*1e3)
            data_trig_ieeg_raw = np.array([data_raw[f"{chan}x00"][0] for chan in chan_trig_extraction]).reshape(-1)/data_trig_ieeg_srate
            data_trig_ieeg_raw = np.sort(data_trig_ieeg_raw)
            data_trig_ieeg = data_trig_ieeg_raw - trial_start

        else:

            trig_ieeg = data_raw[f"{chan_trig_extraction}x00"].reshape(-1)
            trig_ieeg_srate = int(data_raw[srate_trig_extraction][0,0]*1e3)

            trig_ieeg_sig_thresh = (np.median(trig_ieeg) - np.min(trig_ieeg))*0.25
            trig_ieeg_clean = (trig_ieeg < trig_ieeg_sig_thresh)*-1

            trig_ieeg_clean_i = np.where(np.diff(trig_ieeg_clean) != 0)[0]
            data_trig_ieeg_raw = trig_ieeg_clean_i/trig_ieeg_srate

            diff_trig_i = np.diff(data_trig_ieeg_raw)
            diff_trig_i_log = np.log(diff_trig_i)*-1
            thresh_trial_start_stop = (diff_trig_i_log).max()*thresh_trial_start_stop_allpatient[sujet][f"trial{trial_i}"][0]
            _trig_i_list_fast = np.where(thresh_trial_start_stop < diff_trig_i_log)[0]
            start_ieeg_trig = data_trig_ieeg_raw[_trig_i_list_fast[np.argmax(np.diff(_trig_i_list_fast))]+1]
            stop_ieeg_trig = data_trig_ieeg_raw[-1]
            
            if debug:

                plt.plot(trig_ieeg)
                plt.show()

                plt.plot(diff_trig_i_log)
                plt.hlines(thresh_trial_start_stop, xmin=0, xmax=diff_trig_i_log.size, color='r')
                plt.show()

                plt.plot(trig_ieeg_clean)
                plt.vlines(trig_ieeg_clean_i, ymin=trig_ieeg_clean.min(), ymax=trig_ieeg_clean.max(), color='r')
                plt.vlines([start_ieeg_trig*trig_ieeg_srate, stop_ieeg_trig*trig_ieeg_srate], ymin=trig_ieeg_clean.min()*1.2, ymax=trig_ieeg_clean.max()*1.2, color='g')
                plt.show()

        data_trig_ieeg_diff = np.diff(data_trig_ieeg_raw, append=data_trig_ieeg_raw[-1]+10)

        clean_ieeg_thresh_sec = 0.2
        sel_vec_clean_short = data_trig_ieeg_diff < clean_ieeg_thresh_sec 
        data_trig_ieeg_cleaned_short = data_trig_ieeg_raw[~sel_vec_clean_short]
        trig_to_chunk_ieeg = (np.array([data_trig_ieeg_cleaned_short[0], data_trig_ieeg_cleaned_short[-1]])*srate_ieeg).astype(int)

        data_ieeg_chunked = data_ieeg[:,trig_to_chunk_ieeg[0]:trig_to_chunk_ieeg[-1]]
        data_trig_ieeg = data_trig_ieeg_cleaned_short - data_trig_ieeg_cleaned_short[0]

        if debug:

            time_vec_ieeg = np.arange(data_ieeg.shape[-1])/srate_ieeg
            sig = data_ieeg[0]

            plt.plot(time_vec_ieeg, sig)
            plt.vlines(data_trig_ieeg_cleaned_short, ymin=sig.min(), ymax=sig.max(), color='r', label='raw')
            plt.vlines(data_trig_ieeg_cleaned_short[sel_vec_clean_short], ymin=sig.min()*-2, ymax=sig.max()*2, color='g', label='excluded')
            plt.legend()
            plt.show()

            time_vec_ieeg_chunk = np.arange(data_ieeg_chunked.shape[-1])/srate_ieeg
            sig = data_ieeg_chunked[0]

            plt.plot(time_vec_ieeg_chunk, sig)
            plt.vlines(data_trig_ieeg, ymin=sig.min(), ymax=sig.max(), color='r', label='raw')
            plt.legend()
            plt.show()

        #### SL
        resp_clean = physio.preprocess(resp_raw, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
        resp_chunk = resp_clean[trigger_SL_trial_chunk[trial_i][0]:trigger_SL_trial_chunk[trial_i][-1]]
        trigger_SL_chunk = trigger_SL[(trigger_SL <= trigger_SL_trial_chunk[trial_i][-1]) & (trigger_SL_trial_chunk[trial_i][0] <= trigger_SL)]
        trigger_trial_start_SL = trigger_SL_trial_chunk[trial_i][0]/srate_SL
        trigger_SL_trial = trigger_SL_chunk/srate_SL - trigger_trial_start_SL

        if sujet == 'P04' and trial_i == 1:
            crop_to_correct = trigger_SL_trial[1]
            trigger_SL_trial = trigger_SL_trial[1:] - trigger_SL_trial[1] 
            resp_chunk = resp_chunk[int(crop_to_correct*srate_SL):] 
            trigger_SL_trial[-1] = resp_chunk.size/srate_SL

        if debug:

            time_vec = np.arange(resp_chunk.size)/srate_SL
            plt.plot(time_vec, resp_chunk)
            plt.vlines(trigger_SL_trial, ymin=resp_chunk.min(), ymax=resp_chunk.max(), color='r', label='SL')
            # plt.vlines(data_trig_ieeg, ymin=resp_chunk.min(), ymax=resp_chunk.max(), color='g', label='iEEG')
            plt.legend()
            plt.show()

        #### get the shift
        if debug:

            shift = shift_all_patient[sujet][f"trial{trial_i}"]
            shift = 0

            time_vec = np.arange(resp_chunk.size)/srate_SL
            plt.plot(time_vec, resp_chunk)
            plt.vlines(trigger_SL_trial, ymin=resp_chunk.min()*1.1, ymax=resp_chunk.max()*1.1, color='r', label='SL')
            plt.vlines(data_trig_ieeg+shift, ymin=resp_chunk.min()*0.9, ymax=resp_chunk.max()*0.9, color='g', label='iEEG')
            plt.legend()
            plt.show()

        shift = shift_all_patient[sujet][f"trial{trial_i}"]

        data_trig_ieeg_corrected = data_trig_ieeg+shift
        exclude_before_zero = 0 <= data_trig_ieeg_corrected
        data_trig_ieeg_corrected = data_trig_ieeg_corrected[exclude_before_zero]

        data_trig_ieeg_corrected_final = []
        trigger_SL_trial_corrected_final = []
        window_exclusion = 0.8

        for trig_i, trig_val in enumerate(trigger_SL_trial):

            _win_start, _win_stop = trig_val-window_exclusion/2, trig_val+window_exclusion/2
            ieeg_trig_window_sel = data_trig_ieeg_corrected[(_win_start < data_trig_ieeg_corrected) & (data_trig_ieeg_corrected < _win_stop)]

            if ieeg_trig_window_sel.size != 0:
                _argmin_trig_ieeg = np.argmin(np.abs(ieeg_trig_window_sel - trig_val))
                data_trig_ieeg_corrected_final.append(ieeg_trig_window_sel[_argmin_trig_ieeg])
                trigger_SL_trial_corrected_final.append(trig_val)

            if debug:

                time_vec = np.arange(resp_chunk.size)/srate_SL
                plt.plot(time_vec, resp_chunk)
                plt.vlines(trigger_SL_trial, ymin=resp_chunk.min(), ymax=resp_chunk.max(), color='r', label='SL')
                plt.vlines(data_trig_ieeg_corrected, ymin=resp_chunk.min(), ymax=resp_chunk.max(), color='g', label='iEEG')
                plt.vlines([_win_start, _win_stop], ymin=resp_chunk.min()*-2, ymax=resp_chunk.max()*2, color='b', label='window')
                plt.vlines([trig_val], ymin=resp_chunk.min()*-2, ymax=resp_chunk.max()*2, color='b', label='trig_sel', linestyle='--')
                plt.legend()
                plt.show()

        data_trig_ieeg_corrected_final = np.array(data_trig_ieeg_corrected_final)
        trigger_SL_trial_corrected_final = np.array(trigger_SL_trial_corrected_final)

        if debug:
        
            time_vec = np.arange(resp_chunk.size)/srate_SL
            plt.plot(time_vec, resp_chunk)
            plt.vlines(data_trig_ieeg_corrected_final, ymin=resp_chunk.min(), ymax=resp_chunk.max(), color='r', label='iEEG')
            plt.vlines(trigger_SL_trial_corrected_final, ymin=resp_chunk.min(), ymax=resp_chunk.max(), color='g', label='SL')
            plt.legend()
            plt.show()

        #### remove SL trig if SL longer than iEEG
        trigger_SL_trial_corrected_final = trigger_SL_trial_corrected_final[trigger_SL_trial_corrected_final < data_ieeg.shape[-1]/srate_ieeg]

        if trigger_SL_trial_corrected_final.size != data_trig_ieeg_corrected_final.size:
            raise ValueError('!!! trig not same size)')

        #### resample
        data_ieeg_sync = data_ieeg[:,int(data_trig_ieeg_corrected_final.min()*srate_ieeg):int(data_trig_ieeg_corrected_final.max()*srate_ieeg)]
        data_trig_ieeg_sync = data_trig_ieeg_corrected_final - data_trig_ieeg_corrected_final.min()

        resp_resampled, ieeg_resampled = sync_resp_ieeg(resp_chunk, data_ieeg_sync, trigger_SL_trial_corrected_final, data_trig_ieeg_sync, resp_srate=srate_SL, ieeg_srate=srate_ieeg, target_srate=srate)

        if debug:

            time_vec = np.arange(resp_resampled.shape[-1])/srate
            plt.plot(time_vec, resp_resampled)
            plt.vlines(data_trig_ieeg_sync, ymin=resp_resampled.min(), ymax=resp_resampled.max(), color='r')
            plt.show()

        #### aggregates
        trigger_synced = data_trig_ieeg_sync

        #### process resp
        resp_clean_smooth = physio.smooth_signal(resp_resampled, srate, win_shape='gaussian', sigma_ms=40.0)

        if debug:
            resp_resampled_clean = physio.preprocess(resp_resampled, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
            resp_clean_smooth = physio.smooth_signal(resp_resampled_clean, srate, win_shape='gaussian', sigma_ms=40.0)
            
            plt.plot(resp_resampled, label='raw')
            plt.plot(resp_resampled_clean, label='resp_clean')
            plt.plot(resp_clean_smooth, label='resp_clean_smooth')
            plt.legend()
            plt.show()

        cycles_raw = physio.detect_respiration_cycles(resp_resampled, srate, method="crossing_baseline", baseline_mode='zero')

        #### verify cycle detection
        if debug:
                                
            time_vec = np.arange(resp_resampled.shape[0])/srate
            ymin, ymax = resp_resampled.min(), resp_resampled.max()
            fig, ax = plt.subplots(figsize=(18, 10))
            ax.plot(time_vec, resp_resampled)
            ax.scatter(cycles_raw[:,0]/srate, resp_resampled[cycles_raw[:,0]], color='g', label='inspi_selected')
            ax.scatter(cycles_raw[:,1]/srate, resp_resampled[cycles_raw[:,1]], color='c', label='expi_selected', marker='s')
            ax.scatter(cycles_raw[:,2]/srate, resp_resampled[cycles_raw[:,2]], color='r', label='next_inspi_selected', marker='s')
            ax.vlines(trigger_synced, ymin=ymin, ymax=ymax, color='r')
            plt.legend()
            plt.title(f"{sujet} {trial_i} resp")
            plt.rcParams.update({'font.size': 12})
            plt.show()

        resp_features_trial = physio.compute_respiration_cycle_features(resp_resampled, srate, cycles_raw, baseline=None)

        #### export resp figure
        time_vec = np.arange(resp_resampled.shape[0])/srate

        fig_final_resp, ax = plt.subplots(figsize=(18, 10))
        ax.plot(time_vec, resp_resampled)
        ax.scatter(cycles_raw[:,0]/srate, resp_resampled[cycles_raw[:,0]], color='g', label='inspi_selected')
        ax.scatter(cycles_raw[:,1]/srate, resp_resampled[cycles_raw[:,1]], color='c', label='expi_selected', marker='s')
        plt.legend()
        plt.title(f"{sujet} trial{trial_i+1} resp")
        plt.rcParams.update({'font.size': 12})
        # plt.show()
        plt.close()

        fig_final_resp.savefig(os.path.join(path_export_resp_fig, f"{sujet}_trial{trial_i+1}_resp_detect.png"))

        #### export cycle oc_id for patient
        sujet_oc_id = []
        time_vec_resp = np.arange(resp_resampled.size)/srate
        ctrl_oc_append = 0
        oc_trigger = []

        show_cycle_plot = False

        #row_i, row_val = 13, resp_features_trial.loc[13]
        for row_i, row_val in resp_features_trial.iterrows():

            cycle_start, cycle_stop = row_val['inspi_time'], row_val['next_inspi_time']
            trigger_cycle = trigger_synced[(cycle_start < trigger_synced) & (trigger_synced < cycle_stop)]

            if ctrl_oc_append == 1:

                _oc_id = 1

                if np.sum(row_val['expi_time'] < trigger_cycle) == 2:
                    ctrl_oc_append = 1
                else:
                    ctrl_oc_append = 0

            else:

                _oc_id = 0

                if np.sum(row_val['expi_time'] < trigger_cycle) == 2:
                    _oc_id = 0
                    ctrl_oc_append = 1

                if (trigger_cycle < row_val['expi_time']).any():
                    _oc_id = 2
                    oc_trigger.append(cycle_start)

                if 4 < trigger_cycle.size:
                    raise ValueError('!!! error in cycle trig number !!!')                    

            sujet_oc_id.append(_oc_id)

            if show_cycle_plot:
                resp_chunk_cycle = resp_resampled[int(row_val['inspi_index']):int(row_val['next_inspi_index'])] 
                time_vec_cycle_chunk = time_vec_resp[int(row_val['inspi_index']):int(row_val['next_inspi_index'])] 
                plt.plot(time_vec_cycle_chunk, resp_chunk_cycle)
                plt.vlines(trigger_cycle, ymin=resp_chunk_cycle.min(), ymax=resp_chunk_cycle.max(), color='r')
                plt.title(f"cycle_i : {row_i}, class : {_oc_id}")
                plt.show()

        if len(sujet_oc_id) != resp_features_trial.shape[0]:
            raise ValueError('!!! problem in oc_id size and resp feature size !!!')

        oc_trigger = np.array(oc_trigger)

        #### verify OC trigger
        if debug:

            plt.plot(time_vec_resp, resp_resampled)
            plt.vlines(oc_trigger, ymin=resp_resampled.min(), ymax=resp_resampled.max(), color='r')
            plt.title(f"OC")
            plt.show()

        #### get oc_ratio
        oc_ratio = []
        oc_val = []

        #row_i, row_val = 1, resp_features_trial.loc[1]
        for row_i, row_val in resp_features_trial.iterrows():

            cycle_start, cycle_stop = row_val['inspi_time'], row_val['next_inspi_time']
            resp_cycle_sig = resp_resampled[int(cycle_start*srate):int(cycle_stop*srate)]
            expi_index = int(row_val['expi_index']) - int(row_val['inspi_index'])

            _oc_ratio, oc_top_val = get_oc_values(resp_cycle_sig, expi_index)
            
            oc_ratio.append(_oc_ratio)
            oc_val.append(oc_top_val)

            if debug:

                time_vec = np.arange(resp_cycle_sig.size)[:-1]
                sig_diff = np.diff(resp_cycle_sig)
                expi_point = int(row_val['expi_index']) - int(row_val['inspi_index'])
                mid_inspi = int(expi_point/2)
                start_inspi_dec, start_inspi_asc = int(mid_inspi-expi_point*0.375), int(mid_inspi+expi_point*0.375)
                median_inspi_trough = np.median(resp_cycle_sig[start_inspi_dec:start_inspi_asc])
                oc_top_i = np.argmin(resp_cycle_sig[:expi_point])

                plt.plot(time_vec, resp_cycle_sig[:-1])
                plt.vlines([expi_point], ymin=resp_cycle_sig.min(), ymax=resp_cycle_sig.max(), colors='k')
                plt.vlines([mid_inspi], ymin=resp_cycle_sig.min(), ymax=resp_cycle_sig.max(), colors='k', linestyle='--')
                plt.vlines([start_inspi_dec, start_inspi_asc], ymin=resp_cycle_sig.min(), ymax=resp_cycle_sig.max(), colors='g')
                plt.hlines([median_inspi_trough], xmin=0, xmax=expi_point, colors='b')
                plt.scatter([oc_top_i], resp_cycle_sig[oc_top_i], color='r')
                plt.plot(time_vec,sig_diff)
                plt.show()

        if len(oc_ratio) != resp_features_trial.shape[0]:
            raise ValueError('!!! problem in oc_ratio size and resp feature size !!!')

        #### append resp_features
        resp_features_trial.insert(0, 'sujet', [sujet]*resp_features_trial.shape[0])
        resp_features_trial.insert(1, 'cond', [trial_cond]*resp_features_trial.shape[0])
        resp_features_trial.insert(2, 'trial', [trial_i]*resp_features_trial.shape[0])
        resp_features_trial.insert(3, 'cycle', np.arange(resp_features_trial.shape[0]))
        resp_features_trial.insert(4, 'oc_id', sujet_oc_id)
        resp_features_trial.insert(5, 'oc_ratio', oc_ratio)
        resp_features_trial.insert(6, 'oc_val', oc_val)

        #### trial export
        data_ieeg_alltrial_raw.append(ieeg_resampled)
        resp_alltrial.append(resp_resampled)
        oc_trigger_alltrial.append(oc_trigger)
        respfeature_alltrial.append(resp_features_trial)
        
    #### aggregates
    respfeature_alltrial = pd.concat(respfeature_alltrial)

    return data_ieeg_alltrial_raw, resp_alltrial, oc_trigger_alltrial, respfeature_alltrial, chanlist


def bipolar_ref(data_ieeg_alltrial_raw):

    #### bipolar referencing   
    bipolar_data = []

    #trial_i, trial_cond = 0, trial_list_dict[sujet][0]
    for trial_i, trial_cond in enumerate(trial_list_dict[sujet]):

        _data = data_ieeg_alltrial_raw[trial_i]
                
        # construct vec
        chan_i_list = np.arange(_data.shape[0])
        chan_i_list = chan_i_list[chan_i_list != chan_ref_allsujet[sujet]-1]

        pairs = []

        for chan_i in chan_i_list:

            if chan_i != chan_i_list.size:

                if chan_i != bipolarize_chan_i_lim-1:

                    pairs.append([chan_i, chan_i+1])

        data_bi = [_data[pair[0]] - _data[pair[1]] for pair in pairs]
        bipolar_data.append(np.array(data_bi))
        
    data_ieeg_alltrial_bi = bipolar_data
    chanlist_bi = [f"{_pair[0]}_{_pair[-1]}" for _pair in pairs]

    return data_ieeg_alltrial_bi, chanlist_bi

    



################################
######## LOCAFILE ########
################################

def generate_locafile(sujet, chanlist, filename):

    df_locafile_raw = pd.DataFrame({'sujet' : [sujet]*len(chanlist), 'chan' : chanlist, 'loca_raw' : ['']*len(chanlist), 
                       'SOZ' : np.zeros((len(chanlist))), 'OUT' : np.zeros((len(chanlist))), 'NOISY' : np.zeros((len(chanlist))), 
                       'verified' : np.zeros((len(chanlist))), 'select' : np.ones((len(chanlist))), 'loca_modified' : ['']*len(chanlist),})

    df_locafile_raw.to_excel(os.path.join(path_prep, 'anatomy', sujet, filename))






################################
######## COMPARISON ########
################################


# to compare during preprocessing
#data_pre, data_post, chan_name = raw_eeg.get_data(), raw_preproc_wb_clean.get_data(), 'C3'
def compare_pre_post(data_pre, data_post, chanlist, srate, chan_name):

    # compare before after
    nchan_i = chanlist.index(chan_name)
    x_pre = data_pre[nchan_i,:]
    x_post = data_post[nchan_i,:]
    time = np.arange(x_pre.shape[0]) / srate

    nwind = int(10*srate)
    nfft = nwind
    noverlap = np.round(nwind/2)
    hannw = scipy.signal.windows.hann(nwind)

    hzPxx, Pxx_pre = scipy.signal.welch(x_pre,fs=srate,window=hannw,nperseg=nwind,noverlap=noverlap,nfft=nfft)
    hzPxx, Pxx_post = scipy.signal.welch(x_post,fs=srate,window=hannw,nperseg=nwind,noverlap=noverlap,nfft=nfft)

    plt.plot(time, x_pre, label='x_pre')
    plt.plot(time, x_post, label='x_post')
    plt.title(chan_name)
    plt.legend()

    plt.show()

    plt.semilogy(hzPxx, Pxx_pre, label='Pxx_pre')
    plt.semilogy(hzPxx, Pxx_post, label='Pxx_post')
    plt.title(chan_name)
    plt.legend()
    # plt.xlim(60,360)
    plt.show()









################################
######## PREPROCESSING ########
################################

#data, chan_list, srate, prep_step, trial_i, ref_type = data_alltrial_sel[trial_i], chanlist_sel, srate, prep_step_wb, trial_i, ref_type
def preprocessing_ieeg(data, chan_list, srate, prep_step, trial_i, ref_type):

    print('#### PREPROCESSING ####')

    #### export raw PRE PSD
    fig, ax = plt.subplots(figsize=(10,8))
    for chan_i, chan in enumerate(chan_list):
        x = data[chan_i]
        frex, Pxx = scipy.signal.welch(x, fs=srate, nperseg=int(srate * 4))
        ax.semilogy(frex, Pxx, label=chan)
    plt.legend()
    # plt.show()
    fig.savefig(os.path.join(path_results, 'Pxx', 'raw_preproc', f"{sujet}_PRE_t{trial_i}_{ref_type}.png"))


    #### preproc
    
    # 1. Generate raw structures

    ch_types = ['seeg'] * (np.size(data,0)) # ‘ecg’, ‘stim’, ‘eog’, ‘misc’, ‘seeg’, ‘eeg’

    info = mne.create_info(chan_list, srate, ch_types=ch_types, verbose='critical')
    raw_init = mne.io.RawArray(data, info, verbose='critical')

    if debug:
        for nchan in range(np.size(raw_init.get_data(),0)):

            nwind = int( 10*srate ) 
            nfft = nwind 
            noverlap = np.round(nwind/2) 
            hannw = scipy.signal.windows.hann(nwind)

            hzPxx, Pxx = scipy.signal.welch(raw_init.get_data()[nchan,:], fs=srate, window=hannw, nperseg=nwind, noverlap=noverlap, nfft=nfft)

            plt.semilogy(hzPxx, Pxx)
            plt.title(chan_list[nchan])
            #plt.show()
        #plt.show()
        plt.close()

    del data

    # 2. Initiate preprocessing step



    def mean_centered(raw):
        
        data = raw.get_data()
        
        # mean centered
        data_mc = np.zeros((np.size(data,0),np.size(data,1)))
        for chan in range(np.size(data,0)):
            data_mc[chan,:] = data[chan,:] - np.mean(data[chan,:])
            #### no detrend to keep low derivation
            #data_mc[chan,:] = scipy.signal.detrend(data_mc[chan,:]) 

        # fill raw
        for chan in range(np.size(data,0)):
            raw[chan,:] = data_mc[chan,:]

        del data_mc    

        # verif
        if debug:
            # all data
            duration = .5
            n_chan = 10
            raw.plot(scalings='auto',duration=duration,n_channels=n_chan) # verify

            # compare before after
            compare_pre_post(raw.get_data(), raw.get_data(), 4)


        return raw



    def line_noise_removing(raw):

        if debug:
            raw_post = raw.copy()
        else:
            raw_post = raw

        # raw_post.notch_filter(linenoise_freq, verbose='critical')

        alpha_vec_sel = np.linspace(0,1,100)

        data_notch = raw.get_data()

        #chan_i = 1
        for chan_i in range(data_notch.shape[0]):

            x = data_notch[chan_i]

            frex, Pxx = scipy.signal.welch(x, fs=srate, nperseg=int(srate * 4))

            coef = np.polyfit(frex, np.log(Pxx), deg=5)
            trend = np.polyval(coef, frex)
            Pxx_flat = np.log(Pxx) - trend

            if debug:
                plt.plot(np.log(Pxx))
                plt.plot(trend)
                plt.plot(np.log(Pxx) - trend)
                plt.show()

            peaks, _ = scipy.signal.find_peaks(Pxx_flat, height=np.median(Pxx_flat)+Pxx_flat.std()*3, width=.1, distance=int(frex[-1]/20))
            peaks_filtered_i = peaks[15 < frex[peaks]]
            peaks_filtered = frex[peaks_filtered_i]

            if debug:

                plt.plot(frex, Pxx_flat)
                plt.hlines(np.median(Pxx_flat), xmin=0, xmax=frex[-1], color='r')
                plt.hlines(np.median(Pxx_flat)+Pxx_flat.std()*3, xmin=0, xmax=frex[-1], color='r', linestyle='--')
                plt.scatter(peaks_filtered, Pxx_flat[peaks_filtered_i], color='r')
                plt.show()

            x_notch = x.copy().astype('float')
            x_full_notch = x.copy().astype('float')

            for freq_notch_i, freq_notch in zip(peaks_filtered_i, peaks_filtered):

                b, a = scipy.signal.iirnotch(w0=freq_notch, Q=100, fs=srate)

                sos = scipy.signal.tf2sos(b, a)
                freq_x_full_notch = scipy.signal.sosfiltfilt(sos, x_notch)
                x_full_notch = scipy.signal.sosfiltfilt(sos, x_full_notch)

                best_alpha_extraction = []

                for alpha in alpha_vec_sel:

                    x_notch_test = ((1 - alpha) * x_notch + alpha * freq_x_full_notch)
                    frex, Pxx_test_notch = scipy.signal.welch(x_notch_test, fs=srate, nperseg=int(srate * 4))
                    diff = np.abs(trend[freq_notch_i] - np.log(Pxx_test_notch[freq_notch_i]))
                    best_alpha_extraction.append(diff)

                alpha_sel = alpha_vec_sel[np.argmin(best_alpha_extraction)]
                x_notch = ((1 - alpha_sel) * x_notch + alpha_sel * freq_x_full_notch)

                data_notch[chan_i] = x_notch

            if debug:

                frex, Pxx_full_notch = scipy.signal.welch(x_full_notch, fs=srate, nperseg=int(srate * 4))
                frex, Pxx_notch = scipy.signal.welch(x_notch, fs=srate, nperseg=int(srate * 4))
                plt.semilogy(frex, Pxx_full_notch)
                plt.semilogy(frex, Pxx_notch)
                plt.semilogy(frex, Pxx)
                plt.show()

                plt.plot(x_notch, label='notch')
                plt.plot(x_full_notch, label='full notch')
                plt.plot(x, label='raw')
                plt.legend()
                plt.show()

        # fill raw
        for chan in range(np.size(data_notch,0)):
            raw[chan,:] = data_notch[chan,:]

        if debug:

            # compare before after
            compare_pre_post(raw.get_data(), raw_post.get_data(), chanlist, srate, '1_2')

        return raw_post





    def high_pass(raw, h_freq, l_freq):

        if debug:
            raw_post = raw.copy()
        else:
            raw_post = raw

        #filter_length = int(srate*10) # give sec
        filter_length = 'auto'

        if debug:
            h = mne.filter.create_filter(raw_post.get_data(), srate, l_freq=l_freq, h_freq=h_freq, filter_length=filter_length, method='fir', phase='zero-double', fir_window='hamming', fir_design='firwin2')
            flim = (0.1, srate / 2.)
            mne.viz.plot_filter(h, srate, freq=None, gain=None, title=None, flim=flim, fscale='log')

        raw_eeg_mc_hp = raw_post.filter(l_freq, h_freq, filter_length=filter_length, method='fir', phase='zero-double', fir_window='hamming', fir_design='firwin2', verbose='critical')

        if debug:
            duration = 60.
            n_chan = 20
            raw_eeg_mc_hp.plot(scalings='auto',duration=duration,n_channels=n_chan) # verify

        return raw_post


    



    def low_pass(raw, h_freq, l_freq):

        if debug:
            raw_post = raw.copy()
        else:
            raw_post = raw

        filter_length = int(srate*10) # in samples

        if debug:
            h = mne.filter.create_filter(raw_post.get_data(), srate, l_freq=l_freq, h_freq=h_freq, filter_length=filter_length, method='fir', phase='zero-double', fir_window='hamming', fir_design='firwin2')
            flim = (0.1, srate / 2.)
            mne.viz.plot_filter(h, srate, freq=None, gain=None, title=None, flim=flim, fscale='log')

        raw_post = raw_post.filter(l_freq, h_freq, filter_length=filter_length, method='fir', phase='zero-double', fir_window='hann', fir_design='firwin2', verbose='critical')

        if debug:
            duration = .5
            n_chan = 10
            raw_post.plot(scalings='auto',duration=duration,n_channels=n_chan) # verify


        return raw_post



    def average_reref(raw):

        if debug:
            raw_post = raw.copy()
        else:
            raw_post = raw

        raw_post.set_eeg_reference('average')

        if debug:
            duration = .5
            n_chan = 10
            raw_post.plot(scalings='auto',duration=duration,n_channels=n_chan) # verify


        return raw_post



    # 3. Execute preprocessing 

    if debug:
        raw = raw_init.copy() # first data
    else:
        raw = raw_init

    if prep_step['mean_centered']['execute']:
        print('mean_centered')
        raw_post = mean_centered(raw)
        #compare_pre_post(raw.get_data(), raw_post.get_data(), 5)
        raw = raw_post

    if prep_step['line_noise_removing']['execute']:
        print('line_noise_removing')
        raw_post = line_noise_removing(raw)
        #compare_pre_post(raw.get_data(), raw_post.get_data(), 5)
        raw = raw_post

    if prep_step['high_pass']['execute']:
        print('high_pass')
        h_freq = prep_step['high_pass']['params']['h_freq']
        l_freq = prep_step['high_pass']['params']['l_freq']
        raw_post = high_pass(raw, h_freq, l_freq)
        #compare_pre_post(raw.get_data(), raw_post.get_data(), 5)
        raw = raw_post

    if prep_step['low_pass']['execute']:
        print('low_pass')
        h_freq = prep_step['low_pass']['params']['h_freq']
        l_freq = prep_step['low_pass']['params']['l_freq']
        raw_post = low_pass(raw, h_freq, l_freq)
        #compare_pre_post(raw.get_data(), raw_post.get_data(), 5)
        raw = raw_post

    if prep_step['average_reref']['execute']:
        print('average_reref')
        raw_post = average_reref(raw)
        #compare_pre_post(raw.get_data(), raw_post.get_data(), 5)
        raw = raw_post

    data_preproc = raw.get_data()

    del raw, raw_init, raw_post

    #### export raw PSD
    fig, ax = plt.subplots(figsize=(10,8))
    for chan_i, chan in enumerate(chan_list):
        x = data_preproc[chan_i]
        frex, Pxx = scipy.signal.welch(x, fs=srate, nperseg=int(srate * 4))
        ax.semilogy(frex, Pxx, label=chan)
    plt.legend()
    # plt.show()
    fig.savefig(os.path.join(path_results, 'Pxx', 'raw_preproc', f"{sujet}_POST_t{trial_i}_{ref_type}.png"))

    plt.close('all')

    return data_preproc








################################
######## EXPORT ########
################################

def export_data(sujet, data_preproc_alltrial, resp_alltrial, oc_trigger_alltrial, respfeature_alltrial, chanlist, ref_type):

    print('#### SAVE ####')

    #### config
    path_export = os.path.join(path_prep, 'export', sujet)
    path_export_respfeatures = os.path.join(path_precompute, 'respi', 'respfeatures')

    for trial_i, trial in enumerate(trial_list_dict[sujet]):
    
        #### agregate data
        data_all = np.concat([data_preproc_alltrial[trial_i], resp_alltrial[trial_i].reshape(1,-1)], axis=0)
        chanlist_all = chanlist.copy()
        chanlist_all.append('resp')
        xr_dict = {'chan' : chanlist_all, 'time' : np.arange(data_all.shape[-1])/srate}
        xr_data_export = xr.DataArray(data=data_all, dims=xr_dict.keys(), coords=xr_dict)

        #### save
        df_oc_trig = pd.DataFrame({'oc_trig' : oc_trigger_alltrial[trial_i]})

        xr_data_export.to_netcdf(os.path.join(path_export, f"{sujet}_trial{trial_i+1}_preproc_{ref_type}.nc"))    
        df_oc_trig.to_excel(os.path.join(path_export, f"{sujet}_trial{trial_i+1}_oc_trig_{ref_type}.xlsx"))
        

    #### save chanlist + respfeatures
    df_chanlist = pd.DataFrame({'chanlist' : chanlist})

    df_chanlist.to_excel(os.path.join(path_export, f"{sujet}_chanlist_{ref_type}.xlsx"))

    respfeature_alltrial.to_excel(os.path.join(path_export_respfeatures, f"{sujet}_respfeatures.xlsx"))



        




























################################
######## EXECUTE ########
################################


if __name__== '__main__':


    #sujet = sujet_list[0]
    for sujet in sujet_list:

        print(f"{sujet}")

        #ref_type = 'bi'
        for ref_type in ref_type_list:
        
            #### verify computing
            path_export_prep = os.path.join(path_prep, 'export', sujet)

            if os.path.exists(os.path.join(path_export_prep, f'{sujet}_trial1_preproc_{ref_type}.nc')):
                print('ALREADY COMPUTED')
                continue


            ################################
            ######## EXTRACT DATA ########
            ################################

            data_ieeg_alltrial_raw, resp_alltrial, oc_trigger_alltrial, respfeature_alltrial, chanlist = extract_data(sujet)   

            if ref_type == 'bi':
                data_ieeg_alltrial_raw, chanlist = bipolar_ref(data_ieeg_alltrial_raw)

            #### export
            path_export = os.path.join(path_prep, 'export', sujet)

            for trial_i, trial in enumerate(trial_list_dict[sujet]):
            
                #### agregate data
                data_all = np.concat([data_ieeg_alltrial_raw[trial_i], resp_alltrial[trial_i].reshape(1,-1)], axis=0)
                chanlist_all = chanlist.copy()
                chanlist_all.append('resp')
                xr_dict = {'chan' : chanlist_all, 'time' : np.arange(data_all.shape[-1])/srate}
                xr_data_export = xr.DataArray(data=data_all, dims=xr_dict.keys(), coords=xr_dict)
        
                #### save
                xr_data_export.to_netcdf(os.path.join(path_export, f"{sujet}_trial{trial_i+1}_raw_{ref_type}.nc"))


            ################################
            ######## LOCAFILE ########
            ################################

            generate_locafile(sujet, chanlist, f"{sujet}_locafile_raw_{ref_type}.xlsx")
        

            ########################################
            ######## PREPROCESSING & EXPORT ########
            ########################################

            data_preproc_alltrial = []

            for trial_i, trial in enumerate(trial_list_dict[sujet]):

                data_preproc = preprocessing_ieeg(data_ieeg_alltrial_raw[trial_i], chanlist, srate, prep_step_wb, trial_i, ref_type)
                data_preproc_alltrial.append(data_preproc)

            export_data(sujet, data_preproc_alltrial, resp_alltrial, oc_trigger_alltrial, respfeature_alltrial, chanlist, ref_type)














