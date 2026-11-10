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
def epoch_ERP_wb(sujet, ref_type):

    #### config
    path_export_ERP = os.path.join(path_precompute, 'ERP')
    

    if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_RB_WB_erp_{ref_type}.nc")):
        print(f'{sujet} ALREADY COMPUTED', flush=True)
        return

    #### load
    data_alltrial, oc_list_alltrial = get_alltrials_sujet(sujet, ref_type)
    # respfeatures = get_respfeatures(sujet)
    chanlist, localist = get_chanlocalist(sujet, ref_type)

    respfeatures = get_respfeatures(sujet)

    dict_data_allcond = {cond : [] for cond in conditions}

    #### precompute OC
    #cond = conditions_oc[0]
    for cond in conditions_oc:

        print(f'#### ERP PRECOMPUTE {sujet} {cond} ####', flush=True)

        trial_i_list_cond = [trial_i for trial_i, trial_cond in enumerate(trial_list_dict[sujet]) if trial_cond.find(cond) != -1]

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

            epoch_trial_cond = []

            for epoch_i, epoch_trig in enumerate(oc_list_trial):

                _start, _stop = epoch_trig - chunk_time_prepost_i, epoch_trig + chunk_time_prepost_i

                if _start < 0 or data_ieeg_filt_norm.shape[-1] < _stop:
                    continue

                else:
                    epoch_trial_cond.append(data_ieeg_filt_norm[:,_start:_stop].reshape(chanlist.size, 1, -1))

            epoch_trial_cond = np.concat(epoch_trial_cond, axis=1)
            dict_data_allcond[cond].append(epoch_trial_cond)

            #### RB
            inspi_stats_RB = respfeatures.query(f"cond == '{cond}' and trial == {trial_i} and oc_id == 0")['inspi_index'].values

            epoch_trial_RB = []
            
            for epoch_i, epoch_trig in enumerate(inspi_stats_RB):

                _start, _stop = epoch_trig - chunk_time_prepost_i, epoch_trig + chunk_time_prepost_i

                if _start < 0 or data_ieeg_filt_norm.shape[-1] < _stop:
                    continue

                else:
                    epoch_trial_RB.append(data_ieeg_filt_norm[:,_start:_stop].reshape(chanlist.size, 1, -1))

            epoch_trial_RB = np.concat(epoch_trial_RB, axis=1)

            dict_data_allcond['RB'].append(epoch_trial_RB)
            
    #### save
    for cond in conditions:

        time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
        data_concat = np.concat(dict_data_allcond[cond], axis=1)
        dict_xr_erp_coords = {'chan' : chanlist, 'loca' : ('chan', localist['loca_modified'].values), 'epoch' : np.arange(data_concat.shape[1]), 'time' : time_vec_epoch} 
        xr_erp_export = xr.DataArray(data_concat, dims=("chan", "epoch", "time"), coords=dict_xr_erp_coords)

        xr_erp_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_{cond}_WB_erp_{ref_type}.nc"))

    





#sujet = sujet_list[0]
def epoch_ERP_HFA(sujet, ref_type, norm_token=None):

    #### config
    path_export_ERP = os.path.join(path_precompute, 'ERP')

    if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_RB_HFA_erp_{ref_type}.nc")):
        print(f'{sujet} ALREADY COMPUTED', flush=True)
        return

    #### load
    data_alltrial, oc_list_alltrial = get_alltrials_sujet(sujet, ref_type)
    # respfeatures = get_respfeatures(sujet)
    chanlist, localist = get_chanlocalist(sujet, ref_type)

    respfeatures = get_respfeatures(sujet)

    dict_data_allcond = {cond : [] for cond in conditions}

    #### precompute OC
    #cond = conditions_oc[0]
    for cond in conditions_oc:

        print(f'#### ERP PRECOMPUTE {sujet} {cond} ####', flush=True)

        trial_i_list_cond = [trial_i for trial_i, trial_cond in enumerate(trial_list_dict[sujet]) if trial_cond.find(cond) != -1]

        #trial_i = trial_i_list_cond[0]
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

            #### rscore
            _med = np.median(data_hfa, axis=1)
            _mad = scipy.stats.median_abs_deviation(data_hfa, axis=1)

            data_ieeg_filt_norm = (data_hfa - _med[:, None]) * 0.6745 / _mad[:, None]

            if debug:

                chan_i = 0
                plt.plot(data_ieeg.data[chan_i], label='raw')
                plt.plot(data_ieeg_filt_norm[chan_i], label='filt_norm')
                plt.legend()
                plt.show()

            #### OC
            oc_list_trial = (oc_list_alltrial[trial_i]['oc_trig'].values*srate).astype('int')
            chunk_time_prepost_i = int(chunk_time_prepost*srate)

            epoch_trial_cond = []

            for epoch_i, epoch_trig in enumerate(oc_list_trial):

                _start, _stop = epoch_trig - chunk_time_prepost_i, epoch_trig + chunk_time_prepost_i

                if _start < 0 or data_ieeg_filt_norm.shape[-1] < _stop:
                    continue

                else:
                    epoch_trial_cond.append(data_ieeg_filt_norm[:,_start:_stop].reshape(chanlist.size, 1, -1))

            epoch_trial_cond = np.concat(epoch_trial_cond, axis=1)
            dict_data_allcond[cond].append(epoch_trial_cond)

            #### RB
            inspi_stats_RB = respfeatures.query(f"cond == '{cond}' and trial == {trial_i} and oc_id == 0")['inspi_index'].values

            epoch_trial_RB = []
            
            for epoch_i, epoch_trig in enumerate(inspi_stats_RB):

                _start, _stop = epoch_trig - chunk_time_prepost_i, epoch_trig + chunk_time_prepost_i

                if _start < 0 or data_ieeg_filt_norm.shape[-1] < _stop:
                    continue

                else:
                    epoch_trial_RB.append(data_ieeg_filt_norm[:,_start:_stop].reshape(chanlist.size, 1, -1))

            epoch_trial_RB = np.concat(epoch_trial_RB, axis=1)

            dict_data_allcond['RB'].append(epoch_trial_RB)
            
    #### save
    for cond in conditions:

        time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
        data_concat = np.concat(dict_data_allcond[cond], axis=1)
        dict_xr_erp_coords = {'chan' : chanlist, 'loca' : ('chan', localist['loca_modified'].values), 'epoch' : np.arange(data_concat.shape[1]), 'time' : time_vec_epoch} 
        xr_erp_export = xr.DataArray(data_concat, dims=("chan", "epoch", "time"), coords=dict_xr_erp_coords)

        xr_erp_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_{cond}_HFA_erp_{ref_type}.nc"))








################################
######## STATS ########
################################


#sujet = sujet_list[0]
def get_perm_stats_wb(sujet, ref_type):

    ######## VERIFY COMPUTATION ########
    if os.path.exists(os.path.join(path_precompute, 'ERP', f"{sujet}_perm_stats_WB_{ref_type}.nc")):
        return

    ######## CONFIG ########
    path_load_ERP = os.path.join(path_precompute, 'ERP')

    chanlist, localist = get_chanlocalist(sujet, ref_type)
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

        baseline_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_RB_WB_erp_{ref_type}.nc"))

        #cond_i, cond = 0, conditions_sel[0]
        for cond_i, cond in enumerate(conditions_sel):

            erp_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_{cond}_WB_erp_{ref_type}.nc"))

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
    xr_cluster_based_perm.to_netcdf(os.path.join(path_precompute, 'ERP', f"{sujet}_perm_stats_WB_{ref_type}.nc"))

    os.chdir(path_memmap)
    try:
        os.remove(f"{sujet}_perm.dat")
        del xr_data
    except:
        pass




#sujet = sujet_list[0]
def get_perm_stats_HFA(sujet, ref_type):

    ######## VERIFY COMPUTATION ########
    if os.path.exists(os.path.join(path_precompute, 'ERP', f"{sujet}_perm_stats_HFA_{ref_type}.nc")):
        return

    ######## CONFIG ########
    path_load_ERP = os.path.join(path_precompute, 'ERP')

    chanlist, localist = get_chanlocalist(sujet, ref_type)
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

        baseline_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_RB_HFA_erp_{ref_type}.nc"))

        #cond_i, cond = 0, conditions_sel[0]
        for cond_i, cond in enumerate(conditions_sel):

            erp_data = xr.open_dataarray(os.path.join(path_load_ERP, f"{sujet}_{cond}_HFA_erp_{ref_type}.nc"))

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
    xr_cluster_based_perm.to_netcdf(os.path.join(path_precompute, 'ERP', f"{sujet}_perm_stats_HFA_{ref_type}.nc"))

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

    #ref_type = 'mn'
    for ref_type in ref_type_list:

        #sujet = sujet_list[3]
        for sujet in sujet_list:

            print(f'#### {sujet} ####', flush=True)
            print(f'COMPUTE', flush=True)

            epoch_ERP_wb(sujet, ref_type)
            epoch_ERP_HFA(sujet, ref_type)
            get_perm_stats_wb(sujet, ref_type)
            get_perm_stats_HFA(sujet, ref_type)

