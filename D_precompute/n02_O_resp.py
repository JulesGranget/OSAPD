
from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import plotly.graph_objects as go
import seaborn as sns
import neo
import h5py

debug = False







########################################
######## COMPUTE STRETCH ########
########################################

#sujet = sujet_list[0]
def extract_stretch_resp(sujet):

    #### config
    path_export_stretch_respi = os.path.join(path_precompute, 'respi', 'stretch')
    
    #### load respi & CO2
    data_alltrial, oc_list_alltrial = get_alltrials_sujet(sujet)
    respfeatures = get_respfeatures(sujet)
    
    #cond = conditions_oc[4]
    for cond in conditions_oc:

        print(f"computing {cond}...")

        trial_i_list_cond = [trial_i for trial_i, trial_cond in enumerate(trial_list_dict[sujet]) if trial_cond.find(cond) != -1]
        alltrial_stretch_resp = []

        #trial_i = trial_i_list_cond[0]
        for trial_i in trial_i_list_cond:

            xr_data = data_alltrial[trial_i]
            resp = xr_data.sel(chan='resp').data

            respfeatures_trial = respfeatures.query(f"trial == {trial_i}")

            if debug:

                plt.plot(scipy.stats.zscore(resp), label='resp')
                plt.legend()
                plt.show()
            
            #### stretch
            resp_stretch, _ = stretch_data(respfeatures_trial, stretch_point_TF, resp, srate)

            if debug:
                for cycle_i in range(resp_stretch.shape[0]):
                    plt.plot(resp_stretch[cycle_i])
                plt.show()

                plt.plot(resp_stretch.mean(axis=0))
                plt.show()

            #### append all trial
            alltrial_stretch_resp.append(resp_stretch)

        #### concat alltrials
        alltrial_stretch_resp = np.concat(alltrial_stretch_resp)

        #### save
        xr_dict_stretch_resp = {'cycle' : np.arange(alltrial_stretch_resp.shape[0]), 'phase' : np.arange(stretch_point_TF)}
        xr_stretch_resp_cond = xr.DataArray(alltrial_stretch_resp, dims=xr_dict_stretch_resp.keys(), coords=xr_dict_stretch_resp)

        xr_stretch_resp_cond.to_netcdf(os.path.join(path_export_stretch_respi, f"{sujet}_{cond}.nc"))









########################################
######## COMPUTE ERP ########
########################################

#sujet = sujet_list[0]
def extract_ERP_resp(sujet):

    #### config
    path_export_ERP = os.path.join(path_precompute, 'respi', 'ERP')

    #### load
    data_alltrial, oc_list_alltrial = get_alltrials_sujet(sujet)

    #### precompute
    #cond = conditions_oc[0]
    for cond in conditions_oc:

        if os.path.exists(os.path.join(path_export_ERP, f"{sujet}_{cond}_ERP_resp.nc")):
            print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
            continue
        else:
            print(f'#### TF PRECOMPUTE {sujet} {cond} ####', flush=True)

        trial_i_list_cond = [trial_i for trial_i, trial_cond in enumerate(trial_list_dict[sujet]) if trial_cond.find(cond) != -1]
        erp_cond = []

        #trial_i = trial_i_list_cond[0]
        for trial_i in trial_i_list_cond:

            #### data
            data_ieeg = data_alltrial[trial_i].sel(chan=data_alltrial[trial_i]['chan'] == 'resp')
            time_vec_trial = data_ieeg['time'].values

            #### rscore
            _med = np.median(data_ieeg.values, axis=1)
            _mad = scipy.stats.median_abs_deviation(data_ieeg.values, axis=1)

            data_ieeg_norm = (data_ieeg - _med[:, None]) * 0.6745 / _mad[:, None]
            data_ieeg_norm = data_ieeg_norm.data.reshape(-1)

            #### respfeatures
            oc_list_trial = oc_list_alltrial[trial_i]['oc_trig'].values

            epoch_trial = []

            for epoch_i, epoch_trig in enumerate(oc_list_trial):

                _start, _stop = int(epoch_trig*srate) - int(chunk_time_prepost*srate), int(epoch_trig*srate) + int(chunk_time_prepost*srate)

                if _start < 0 or time_vec_trial.shape[0] < _stop:
                    continue

                epoch_trial.append(data_ieeg_norm[_start:_stop])

            epoch_trial = np.stack(epoch_trial)
            erp_cond.append(epoch_trial)
            
        #### concat
        erp_cond = np.concat(erp_cond, axis=0)

        #### save stretch tf
        print(f'SAVE ERP {sujet} {cond}', flush=True)

        time_vec_epoch = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
        dict_xr_erp_coords = {'epoch' : np.arange(erp_cond.shape[0]), 'time' : time_vec_epoch} 
        xr_erp_export = xr.DataArray(erp_cond, dims=("epoch", "time"), coords=dict_xr_erp_coords)

        xr_erp_export.to_netcdf(os.path.join(path_export_ERP, f"{sujet}_{cond}_ERP_resp.nc"))







########################################
######## EXPORT MAIN FIGURE ########
########################################

#sujet = sujet_list[0]
def extract_mean_figures_erp(sujet):

    #### config
    path_export_mean_fig = os.path.join(path_results, 'respi', 'mean')
    path_load_stretch_respi = os.path.join(path_precompute, 'respi', 'ERP')
    
    #### load
    erp_resp = {cond : xr.load_dataarray(os.path.join(path_load_stretch_respi, f"{sujet}_{cond}_ERP_resp.nc")) for cond in conditions_oc}
    respfeatures = get_respfeatures(sujet)

    #### export mean fig
    conditions_cycle_plot = ['O']
    time_vec_plot = erp_resp[conditions_oc[0]]['time'].data

    fig_mean, axs = plt.subplots(nrows=2, ncols=len(conditions_cycle_plot), figsize=(10,8))
    
    for col_i, cond_cycle in enumerate(conditions_cycle_plot):

        for row_i, plot_type in enumerate(['mean', 'allcycle']):

            ax = axs[row_i]
            cycle_number_oc = {}

            for cond_oc in conditions_oc:

                _cycle_plot  = erp_resp[cond_oc].data

                cycle_number_oc[cond_oc] = _cycle_plot.shape[0]

                mins, maxs = [], []

                if plot_type == 'mean':
                
                    _mean = _cycle_plot.mean(axis=0)
                    _std = _cycle_plot.std(axis=0)
                    (line,) = ax.plot(time_vec_plot, _mean, label=cond_oc)
                    ax.fill_between(time_vec_plot, _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)
                    min, max = (_mean-_std).min(), (_mean+_std).max()
                    mins.append(min) 
                    maxs.append(max)
                
                elif plot_type == 'allcycle':

                    for cycle_i in range(_cycle_plot.shape[0]):
                        ax.plot(time_vec_plot, _cycle_plot[cycle_i])
                    min, max = _cycle_plot.min(), _cycle_plot.max()
                    mins.append(min) 
                    maxs.append(max)

            ax.legend()
            ax.set_title(f"{cond_cycle} {plot_type} oc300:{cycle_number_oc['oc300']}, oc600:{cycle_number_oc['oc600']}")
            ax.vlines(0, ymin=np.array(mins).min(), ymax=np.array(maxs).max(), color='r')         

    plt.suptitle(f'{sujet} ERP')
    plt.tight_layout()
    plt.rcParams.update({'font.size': 12})
    plt.legend()
    # plt.show()

    fig_mean.savefig(os.path.join(path_export_mean_fig, f"{sujet}_erp_mean.png"))




#sujet = sujet_list[0]
def extract_mean_figures_stretch(sujet):

    #### config
    path_export_mean_fig = os.path.join(path_results, 'respi', 'mean')
    path_load_stretch_respi = os.path.join(path_precompute, 'respi', 'stretch')
    
    #### load
    stretch_resp = {cond : xr.load_dataarray(os.path.join(path_load_stretch_respi, f"{sujet}_{cond}.nc")) for cond in conditions_oc}
    respfeatures = get_respfeatures(sujet)

    #### export mean fig
    fig_mean, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,8))

    oc_id_sel = {'RB' : [0], 'O' : [2]}
    
    for col_i, cond_cycle in enumerate(conditions_cycle):

        for row_i, plot_type in enumerate(['mean', 'allcycle']):

            ax = axs[row_i,col_i]
            cycle_number_oc = {}

            for cond_oc in conditions_oc:

                respfeatures_to_plot = respfeatures.query(f"cond == '{cond_oc}'")
                respfeatures_to_plot['cycle_i_cond'] = np.arange(respfeatures_to_plot.shape[0])
                respfeatures_to_plot = respfeatures_to_plot.query(f"oc_id in {oc_id_sel[cond_cycle]}")
                stretch_cycle_plot  = stretch_resp[cond_oc].data[respfeatures_to_plot['cycle_i_cond'].values]

                cycle_number_oc[cond_oc] = stretch_cycle_plot.shape[0]

                mins, maxs = [], []

                if plot_type == 'mean':
                
                    _mean = stretch_cycle_plot.mean(axis=0)
                    _std = stretch_cycle_plot.std(axis=0)
                    (line,) = ax.plot(_mean, label=cond_oc)
                    ax.fill_between(range(_mean.size), _mean - _std, _mean + _std, color=line.get_color(), alpha=0.1)
                    min, max = (_mean-_std).min(), (_mean+_std).max()
                    mins.append(min) 
                    maxs.append(max)
                
                elif plot_type == 'allcycle':

                    for cycle_i in range(stretch_cycle_plot.shape[0]):
                        ax.plot(stretch_cycle_plot[cycle_i])
                    min, max = stretch_cycle_plot.min(), stretch_cycle_plot.max()
                    mins.append(min) 
                    maxs.append(max)

            ax.legend()
            ax.set_title(f"{cond_cycle} {plot_type} oc300:{cycle_number_oc['oc300']}, oc600:{cycle_number_oc['oc600']}")
            ax.vlines(int(stretch_point_TF/2), ymin=np.array(mins).min(), ymax=np.array(maxs).max(), color='r')         

    plt.suptitle(f'{sujet} stretch')
    plt.tight_layout()
    plt.rcParams.update({'font.size': 12})
    plt.legend()
    # plt.show()

    fig_mean.savefig(os.path.join(path_export_mean_fig, f"{sujet}_stretch_mean.png"))








################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        extract_stretch_resp(sujet)
        extract_ERP_resp(sujet)

        extract_mean_figures_stretch(sujet)
        extract_mean_figures_erp(sujet)
















