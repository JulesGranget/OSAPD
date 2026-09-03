import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import pandas as pd
import joblib
from sklearn import linear_model
import gc

from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *


debug = False





################################
######## EPOCH ERP ########
################################

#sujet = sujet_list[0]
def epoch_ERP_wb(sujet):

    #### config
    path_export_ERP = os.path.join(path_precompute, 'ERP', sujet)

    #### load
    data_alltrial, oc_list_alltrial = get_alltrials_sujet(sujet)
    # respfeatures = get_respfeatures(sujet)
    chanlist, localist = get_chanlocalist(sujet)

    #### precompute OC
    #cond = conditions_oc[0]
    for cond in conditions_oc:

        if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_{cond}_WB_erp.nc")):
            print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
            continue
        else:
            print(f'#### TF PRECOMPUTE {sujet} {cond} ####', flush=True)

        trial_i_list_cond = [trial_i for trial_i, trial_cond in enumerate(trial_list_dict[sujet]) if trial_cond.find(cond) != -1]
        erp_cond = []

        #trial_i = trial_i_list_cond[0]
        for trial_i in trial_i_list_cond:

            #### data
            data_ieeg = data_alltrial[trial_i].sel(chan=data_alltrial[trial_i]['chan'] != 'resp')
            time_vec_trial = data_ieeg['time'].values

            #### detrend + low pass 45Hz
            data_ieeg_filt = []

            for chan_i, chan in enumerate(chanlist):
                x = scipy.signal.detrend(data_ieeg.sel(chan=chan).data, type='linear')
                x = iirfilt(x, srate, lowcut=0.05, highcut=None)
                x = iirfilt(x, srate, lowcut=None, highcut=45)
                data_ieeg_filt.append(x)

            data_ieeg_filt = np.array(data_ieeg_filt)

            #### rscore
            _med = np.median(data_ieeg_filt, axis=1)
            _mad = scipy.stats.median_abs_deviation(data_ieeg_filt, axis=1)

            data_ieeg_filt_norm = (data_ieeg_filt - _med[:, None]) * 0.6745 / _mad[:, None]

            if debug:

                chan_i = 0
                plt.plot(data_ieeg.data[chan_i], label='raw')
                plt.plot(data_ieeg_filt_norm[chan_i], label='filt_norm')
                plt.legend()
                plt.show()

            #### OC
            oc_list_trial = (oc_list_alltrial[trial_i]['oc_trig'].values*srate).astype('int')
            chunk_time_prepost_i = int(chunk_time_prepost*srate)

            epoch_trial = []

            for epoch_i, epoch_trig in enumerate(oc_list_trial):

                _start, _stop = epoch_trig - chunk_time_prepost_i, epoch_trig + chunk_time_prepost_i

                if _start < 0 or data_ieeg_filt_norm.shape[-1] < _stop:
                    continue

                else:
                    epoch_trial.append(data_ieeg_filt_norm[:,_start:_stop].reshape(chanlist.size, 1, -1))

            epoch_trial = np.concat(epoch_trial, axis=1)
            erp_cond.append(epoch_trial)
            
        #### concat
        erp_cond = np.concat(erp_cond, axis=1)

        #### save stretch tf
        print(f'SAVE ERP {sujet} {cond}', flush=True)

        time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
        dict_xr_erp_coords = {'chan' : chanlist, 'loca' : ('chan', localist['loca_modified'].values), 'epoch' : np.arange(erp_cond.shape[1]), 'time' : time_vec_epoch} 
        xr_erp_export = xr.DataArray(erp_cond, dims=("chan", "epoch", "time"), coords=dict_xr_erp_coords)

        xr_erp_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_{cond}_WB_erp.nc"))

    #### precompute RB
    cond = 'RB'
    if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_RB_WB_erp.nc")):
        print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
        return
    else:
        print(f'#### TF PRECOMPUTE {sujet} {cond} ####', flush=True)

    erp_cond = []

    respfeatures = get_respfeatures(sujet)

    #trial_i, trial = 0, trial_list_dict[sujet][0]
    for trial_i, trial in enumerate(trial_list_dict[sujet]):

        #### data
        data_ieeg = data_alltrial[trial_i].sel(chan=data_alltrial[trial_i]['chan'] != 'resp')
        time_vec_trial = data_ieeg['time'].values

        #### detrend + low pass 45Hz
        data_ieeg_filt = []

        for chan_i, chan in enumerate(chanlist):
            x = scipy.signal.detrend(data_ieeg.sel(chan=chan).data, type='linear')
            x = iirfilt(x, srate, lowcut=0.05, highcut=None)
            x = iirfilt(x, srate, lowcut=None, highcut=45)
            data_ieeg_filt.append(x)

        data_ieeg_filt = np.array(data_ieeg_filt)

        #### rscore
        _med = np.median(data_ieeg_filt, axis=1)
        _mad = scipy.stats.median_abs_deviation(data_ieeg_filt, axis=1)

        data_ieeg_filt_norm = (data_ieeg_filt - _med[:, None]) * 0.6745 / _mad[:, None]

        if debug:

            chan_i = 0
            plt.plot(data_ieeg.data[chan_i], label='raw')
            plt.plot(data_ieeg_filt_norm[chan_i], label='filt_norm')
            plt.legend()
            plt.show()

        #### inspi
        inspi_starts_i_list = respfeatures.query(f"trial == {trial_i} and oc_id == 0")['inspi_index'].values
        chunk_time_prepost_i = int(chunk_time_prepost*srate)

        epoch_trial = []

        for epoch_i, epoch_trig in enumerate(inspi_starts_i_list):

            _start, _stop = epoch_trig - chunk_time_prepost_i, epoch_trig + chunk_time_prepost_i

            if _start < 0 or data_ieeg_filt_norm.shape[-1] < _stop:
                continue

            else:
                epoch_trial.append(data_ieeg_filt_norm[:,_start:_stop].reshape(chanlist.size, 1, -1))

        epoch_trial = np.concat(epoch_trial, axis=1)
        erp_cond.append(epoch_trial)
        
    #### concat
    erp_cond = np.concat(erp_cond, axis=1)

    #### save stretch tf
    print(f'SAVE ERP {sujet} {cond}', flush=True)

    time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
    dict_xr_erp_coords = {'chan' : chanlist, 'loca' : ('chan', localist['loca_modified'].values), 'epoch' : np.arange(erp_cond.shape[1]), 'time' : time_vec_epoch} 
    xr_erp_export = xr.DataArray(erp_cond, dims=("chan", "epoch", "time"), coords=dict_xr_erp_coords)

    xr_erp_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_RB_WB_erp.nc"))





#sujet = sujet_list[0]
def epoch_ERP_HFA(sujet):

    #### config
    path_export_ERP = os.path.join(path_precompute, 'ERP', sujet)

    #### load
    data_alltrial, oc_list_alltrial = get_alltrials_sujet(sujet)
    chanlist, localist = get_chanlocalist(sujet)

    for cond in conditions_oc:

        if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_{cond}_HFA_erp.nc")):
            print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
            continue
        else:
            print(f'#### HFA PRECOMPUTE {sujet} {cond} ####', flush=True)

        trial_i_list_cond = [trial_i for trial_i, trial_cond in enumerate(trial_list_dict[sujet]) if trial_cond.find(cond) != -1]

        hfa_cond = []

        for trial_i in trial_i_list_cond:

            #### data
            data_ieeg = data_alltrial[trial_i].sel(chan=data_alltrial[trial_i]['chan'] != 'resp')

            #### HFA
            data_hfa = []

            for chan_i, chan in enumerate(chanlist):

                x = data_ieeg.sel(chan=chan).data

                #### detrend
                x = scipy.signal.detrend(x, type='linear')

                #### bandpass HFA
                x = iirfilt(x, srate, lowcut=HFA_freq[0], highcut=HFA_freq[-1])

                #### Hilbert analytic signal
                x_hilbert = scipy.signal.hilbert(x)

                #### HFA power
                x_power = np.abs(x_hilbert) ** 2

                data_hfa.append(x_power)

            data_hfa = np.array(data_hfa)

            if debug:

                plt.plot(data_hfa[0])
                plt.show()

            #### rscore
            _med = np.median(data_hfa, axis=1)

            _mad = scipy.stats.median_abs_deviation(data_hfa, axis=1)

            _mad[_mad == 0] = np.nan

            data_hfa_norm = ((data_hfa - _med[:, None]) * 0.6745 / _mad[:, None])

            if debug:

                chan_i = 0

                plt.figure(figsize=(12, 5))

                plt.plot(data_hfa[chan_i], label='HFA power')
                plt.plot(data_hfa_norm[chan_i], label='HFA power rscore')
                plt.legend()
                plt.show()

            #### epoch
            oc_list_trial = (oc_list_alltrial[trial_i]['oc_trig'].values * srate).astype(int)

            chunk_time_prepost_i = int(chunk_time_prepost * srate)

            epoch_trial = []

            for epoch_i, epoch_trig in enumerate(oc_list_trial):

                _start = epoch_trig - chunk_time_prepost_i
                _stop = epoch_trig + chunk_time_prepost_i

                if (_start < 0 or data_hfa_norm.shape[-1] < _stop):
                    continue

                epoch_trial.append(data_hfa_norm[:, _start:_stop].reshape(chanlist.size, 1, -1))

            if len(epoch_trial) == 0:
                continue

            epoch_trial = np.concatenate(epoch_trial, axis=1)

            hfa_cond.append(epoch_trial)

        hfa_cond = np.concatenate(hfa_cond, axis=1)

        #### save
        print(f'SAVE HFA {sujet} {cond}', flush=True)

        time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1 / srate)

        dict_xr_hfa_coords = {'chan': chanlist, 'loca': ('chan', localist['loca_modified'].values), 'epoch': np.arange(hfa_cond.shape[1]), 'time': time_vec_epoch}

        xr_hfa_export = xr.DataArray(hfa_cond, dims=("chan", "epoch", "time"), coords=dict_xr_hfa_coords)

        xr_hfa_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_{cond}_HFA_erp.nc"))

    #### RB
    cond = 'RB'

    if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_RB_HFA_erp.nc")):
        print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
        return

    else:

        print(f'#### HFA PRECOMPUTE {sujet} {cond} ####', flush=True)

    hfa_cond = []

    respfeatures = get_respfeatures(sujet)

    for trial_i, trial in enumerate(trial_list_dict[sujet]):

        #### data
        data_ieeg = data_alltrial[trial_i].sel(chan=data_alltrial[trial_i]['chan'] != 'resp')

        time_vec_trial = data_ieeg['time'].values

        #### HFA
        data_hfa = []

        for chan_i, chan in enumerate(chanlist):

            x = data_ieeg.sel(chan=chan).data

            #### detrend
            x = scipy.signal.detrend(x, type='linear')

            #### bandpass HFA
            x = iirfilt(x, srate, lowcut=70, highcut=150)

            #### Hilbert
            x_hilbert = scipy.signal.hilbert(x)

            #### HFA power
            x_power = np.abs(x_hilbert) ** 2

            data_hfa.append(x_power)

        data_hfa = np.array(data_hfa)

        #### rscore

        _med = np.median(data_hfa, axis=1)

        _mad = scipy.stats.median_abs_deviation(data_hfa, axis=1)

        _mad[_mad == 0] = np.nan

        data_hfa_norm = ((data_hfa - _med[:, None]) * 0.6745 / _mad[:, None])

        if debug:

            chan_i = 0

            plt.figure(figsize=(12, 5))

            plt.plot(data_hfa[chan_i], label='HFA power')

            plt.plot(data_hfa_norm[chan_i], label='HFA power rscore')

            plt.legend()
            plt.show()

        #### inspi

        inspi_starts_i_list = respfeatures.query(f"trial == {trial_i} and oc_id == 0")['inspi_index'].values

        chunk_time_prepost_i = int(chunk_time_prepost * srate)

        epoch_trial = []

        for epoch_i, epoch_trig in enumerate(inspi_starts_i_list):

            _start = epoch_trig - chunk_time_prepost_i
            _stop = epoch_trig + chunk_time_prepost_i

            if (_start < 0 or data_hfa_norm.shape[-1] < _stop):
                continue

            epoch_trial.append(data_hfa_norm[:, _start:_stop].reshape(chanlist.size, 1, -1))

        if len(epoch_trial) == 0:
            continue

        epoch_trial = np.concatenate(epoch_trial, axis=1)

        hfa_cond.append(epoch_trial)

    hfa_cond = np.concatenate(hfa_cond, axis=1)

    #### save

    print(f'SAVE HFA {sujet} {cond}', flush=True)

    time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1 / srate)

    dict_xr_hfa_coords = {'chan': chanlist, 'loca': ('chan', localist['loca_modified'].values), 'epoch': np.arange(hfa_cond.shape[1]), 'time': time_vec_epoch}

    xr_hfa_export = xr.DataArray(hfa_cond, dims=("chan", "epoch", "time"), coords=dict_xr_hfa_coords)

    xr_hfa_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_RB_HFA_erp.nc"))








################################
######## STATS ########
################################


#sujet = sujet_list[0]
def get_perm_stats_wb(sujet):

    ######## VERIFY COMPUTATION ########
    if os.path.exists(os.path.join(path_precompute, 'ERP', sujet, f"{sujet}_perm_stats_wb.nc")):
        return

    ######## CONFIG ########
    path_load_ERP = os.path.join(path_precompute, 'ERP', sujet)

    chanlist, localist = get_chanlocalist(sujet)
    conditions_sel = conditions_oc
    time_vec = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)

    ######## INITIATE MEMMAP ########
    path_memmap_computation = os.path.join(path_memmap, f"{sujet}_perm.dat")
    xr_data = np.memmap(path_memmap_computation, dtype=np.float32, mode='w+', shape=(len(conditions_sel), chanlist.size, time_vec.size))

    ######## PARALLELIZATION FUNCTION ########
    #chan_i, chan = 0, chanlist[0]
    def get_cluster_based_perm_one_sujet(chan_i, chan):

        ######## COMPUTE PERMUTATION ########
        print(f'chan{chan} COMPUTE PERMUTATION')

        baseline_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_RB_WB_erp.nc"))

        #cond_i, cond = 0, conditions_sel[0]
        for cond_i, cond in enumerate(conditions_sel):

            erp_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_{cond}_WB_erp.nc"))


            data_baseline = baseline_data.sel(chan=chan).data
            data_cond = erp_data.sel(chan=chan).data

            mask = get_permutation_cluster_1d(data_baseline, data_cond, ERP_n_surrogate, stat_design='between', mode_grouped='median', mode_generate_surr='percentile_time', 
                            mode_select_thresh='percentile_time', percentile_thresh=[2.5, 97.5], size_thresh_alpha=0.01)

            xr_data[cond_i, chan_i, :] = mask

        print(f'{chan} done')

    ######## PARALLELIZATION COMPUTATION ########
    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(get_cluster_based_perm_one_sujet)(chan_i, chan) for chan_i, chan in enumerate(chanlist))

    ######## SAVE ########
    xr_coords = {'cond' : conditions_sel, 'chan' : chanlist, 'time' : time_vec}
    xr_cluster_based_perm = xr.DataArray(data=xr_data, dims=xr_coords.keys(), coords=xr_coords.values())
    xr_cluster_based_perm.to_netcdf(os.path.join(path_precompute, 'ERP', sujet, f"{sujet}_perm_stats_wb.nc"))

    os.chdir(path_memmap)
    try:
        os.remove(f"{sujet}_perm.dat")
        del xr_data
    except:
        pass




#sujet = sujet_list[0]
def get_perm_stats_HFA(sujet):

    ######## VERIFY COMPUTATION ########
    if os.path.exists(os.path.join(path_precompute, 'ERP', sujet, f"{sujet}_perm_stats_HFA.nc")):
        return

    ######## CONFIG ########
    path_load_ERP = os.path.join(path_precompute, 'ERP', sujet)

    chanlist, localist = get_chanlocalist(sujet)
    conditions_sel = conditions_oc
    time_vec = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)

    ######## INITIATE MEMMAP ########
    path_memmap_computation = os.path.join(path_memmap, f"{sujet}_perm.dat")
    xr_data = np.memmap(path_memmap_computation, dtype=np.float32, mode='w+', shape=(len(conditions_sel), chanlist.size, time_vec.size))

    ######## PARALLELIZATION FUNCTION ########
    #chan_i, chan = 0, chanlist[0]
    def get_cluster_based_perm_one_sujet(chan_i, chan):

        ######## COMPUTE PERMUTATION ########
        print(f'chan{chan} COMPUTE PERMUTATION')

        baseline_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_RB_HFA_erp.nc"))

        #cond_i, cond = 0, conditions_sel[0]
        for cond_i, cond in enumerate(conditions_sel):

            erp_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_{cond}_HFA_erp.nc"))

            data_baseline = baseline_data.sel(chan=chan).data
            data_cond = erp_data.sel(chan=chan).data

            mask = get_permutation_cluster_1d(data_baseline, data_cond, ERP_n_surrogate, stat_design='between', mode_grouped='median', mode_generate_surr='percentile_time', 
                            mode_select_thresh='percentile_time', percentile_thresh=[2.5, 97.5], size_thresh_alpha=0.01)

            xr_data[cond_i, chan_i, :] = mask

        print(f'{chan} done')

    ######## PARALLELIZATION COMPUTATION ########
    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(get_cluster_based_perm_one_sujet)(chan_i, chan) for chan_i, chan in enumerate(chanlist))

    ######## SAVE ########
    xr_coords = {'cond' : conditions_sel, 'chan' : chanlist, 'time' : time_vec}
    xr_cluster_based_perm = xr.DataArray(data=xr_data, dims=xr_coords.keys(), coords=xr_coords.values())
    xr_cluster_based_perm.to_netcdf(os.path.join(path_precompute, 'ERP', sujet, f"{sujet}_perm_stats_HFA.nc"))

    os.chdir(path_memmap)
    try:
        os.remove(f"{sujet}_perm.dat")
        del xr_data
    except:
        pass









################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    
    sujet = sujet_list[0]
    for sujet in sujet_list:

        print(f'#### {sujet} ####', flush=True)
        print(f'COMPUTE', flush=True)

        epoch_ERP_wb(sujet)
        epoch_ERP_HFA(sujet)
        get_perm_stats_wb(sujet)
        get_perm_stats_HFA(sujet)

