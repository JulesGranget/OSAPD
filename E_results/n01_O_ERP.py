


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *





################################
######## ERP PLOT ########
################################



def plot_ERP_res(sujet):

    #### config
    path_export_fig = os.path.join(path_results, 'ERP', sujet)

    freq_cond_list = ['WB', 'HFA']
    conditions_sel = ['RB', 'oc300', 'oc600']
    time_vec_ERP = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
    time_vec_plot = (-1 < time_vec_ERP) & (time_vec_ERP < 2)
    path_load = os.path.join(path_precompute, 'ERP', sujet)
    chanlist, localist = get_locachanlist(sujet)

    ### load ERP
    data_ERP = {}

    for freq_type in freq_cond_list:

        data_ERP[freq_type] = {}

        for cond in conditions_sel:        

            data_ERP[freq_type][cond] = xr.open_dataarray(os.path.join(path_load, f"{sujet}_{cond}_{freq_type}_erp.nc"))

    stats_ERP = {}

    for freq_type in freq_cond_list:

        data_ERP[freq_type][cond] = xr.open_dataarray(os.path.join(path_load, f"{sujet}_{cond}_{freq_type}_erp.nc"))

    #### extract vlims
    vlims = {}
    SE = {}

    for row_i, freq_type in enumerate(freq_cond_list):

        vlims[freq_type] = {}
        SE[freq_type] = {}

        for chan_i, chan in enumerate(chanlist):

            _min, _max = [], []
            SE[freq_type][chan] = {}

            for cond in conditions_sel:

                _sig_mat = data_ERP[freq_type][cond].sel(chan=chan).data
                _sig_SE = np.std(_sig_mat, axis=0) / _sig_mat.shape[0] ** 0.5
                _sig = np.median(_sig_mat, axis=0)
                _sig_mSE, _sig_pSE = _sig - _sig_SE, _sig + _sig_SE
                _min.append(_sig_mSE.min())
                _max.append(_sig_pSE.max())

                SE[freq_type][chan][cond] = [_sig_mSE, _sig_pSE]

            vlims[freq_type][chan] = [np.array(_min).min(), np.array(_max).max()]
                
    #### plot
    cond_color = {'RB' : 'g', 'oc300' : 'k', 'oc600' : 'r'}

    #chan_i, chan = 0, chanlist[0]
    for chan_i, chan in enumerate(chanlist):

        fig, axs = plt.subplots(nrows=len(freq_cond_list), ncols=1)

        plt.suptitle(f'{chan}')

        fig.set_figheight(10)
        fig.set_figwidth(10)

        #row_i, freq_type = 0, freq_cond_list[0]
        for row_i, freq_type in enumerate(freq_cond_list):

            ax = axs[row_i]

            #cond = conditions_sel[0]
            for cond in conditions_sel:    

                _sig_mat = data_ERP[freq_type][cond].sel(chan=chan).data
                _sig_med = np.median(_sig_mat, axis=0)
                _sig_mSE, _sig_pSE = SE[freq_type][chan][cond][0], SE[freq_type][chan][cond][1]

                ax.plot(time_vec_ERP[time_vec_plot], _sig_med[time_vec_plot], color=cond_color[cond], label=f"{cond}:c{_sig_mat.shape[0]}")
                ax.fill_between(time_vec_ERP[time_vec_plot], _sig_mSE[time_vec_plot], _sig_pSE[time_vec_plot], color=cond_color[cond], alpha=0.5)

            ax.set_title(f"{freq_type}", fontweight='bold')

            ax.set_ylim(vlims[freq_type][chan][0], vlims[freq_type][chan][1])

            ax.invert_yaxis()

            ax.vlines(0, ymin=vlims[freq_type][chan][0], ymax=vlims[freq_type][chan][1], colors='k')

            ax.legend()

        fig.tight_layout()

        # plt.show()

        #### save
        fig.savefig(os.path.join(path_export_fig, f'{sujet}_{chan}.jpeg'), dpi=150)

        fig.clf()
        plt.close('all')
        gc.collect()

    






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        plot_ERP_res(sujet)

    
