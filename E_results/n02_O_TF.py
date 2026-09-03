



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *






################################
######## TF RELATIVE ########
################################



def tf_relative_patientwise(sujet):

    #### config
    path_load_TF = os.path.join(path_precompute, 'TF', sujet)
    path_export_TF = os.path.join(path_results, 'TF', sujet)
        
    TF_type_list = ['absolute', 'diff']

    #### load extract vlim
    respfeatures = get_respfeatures(sujet)
    chanlist, localist = get_chanlocalist(sujet)

    xr_dict_data = {}
    xr_RB_data = []

    for cond in conditions_oc:

        _xr_load = xr.open_dataarray(os.path.join(path_load_TF, f"{sujet}_{cond}_tf_stretch.nc"))
        oc_id_sel_cond = respfeatures.query(f"cond == '{cond}'").reset_index(drop=True).query(f"oc_id == 2").index.values
        _xr_cond = _xr_load.sel(cycle=oc_id_sel_cond)
        xr_dict_data[cond] = _xr_cond

        oc_id_sel_RB = respfeatures.query(f"cond == '{cond}'").reset_index(drop=True).query(f"oc_id == 0").index.values
        _xr_RB = _xr_load.sel(cycle=oc_id_sel_RB)
        xr_RB_data.append(_xr_RB)

    xr_dict_data['RB'] = xr.concat(xr_RB_data, dim='cycle')

    #### get vlims
    vlims = {}
    
    for TF_type in TF_type_list:

        vlims[TF_type] = {}

        for chan in chanlist:

            vlim = []

            data_baseline = xr_dict_data['RB'].sel(chan=chan).median('cycle')

            for cond in conditions_oc:

                _xr_cond = xr_dict_data[cond].sel(chan=chan).median('cycle')

                if TF_type == 'absolute':
                    
                    vlim.append([_xr_cond.min().values, _xr_cond.max().values])

                elif TF_type == 'relative':

                    _xr_diff = _xr_cond - data_baseline
                    vlim.append([_xr_diff.min().values, _xr_diff.max().values])

            vlim.append([data_baseline.min().values, data_baseline.max().values])
            vlims[TF_type][chan] = np.abs(np.array(vlim)).max()

    #### plot
    phase_vec = xr_dict_data['RB']['phase'].data
    ticks_freq = [2, 8, 12, 30, 50, 60, 150]
    conditions_oc_baseline = ['RB', 'oc300', 'oc600']

    for chan_i, chan in enumerate(chanlist):

        fig, axs = plt.subplots(
            ncols=len(conditions_oc_baseline),
            nrows=len(TF_type_list),
            figsize=(16, 10)
        )

        for row_i, TF_type in enumerate(TF_type_list):

            for col_i, cond in enumerate(conditions_oc_baseline):

                if cond == 'RB':
                    xr_plot = xr_dict_data['RB'].sel(chan=chan).median('cycle').data

                else:
                    if TF_type == 'absolute':
                        xr_plot = xr_dict_data[cond].sel(chan=chan).median('cycle').data

                    elif TF_type == 'diff':
                        xr_plot = xr_dict_data[cond].sel(chan=chan).median('cycle').data - xr_dict_data['RB'].sel(chan=chan).median('cycle').data

                ax = axs[row_i, col_i]

                pcm = ax.pcolormesh(
                    phase_vec,
                    frex,
                    xr_plot,
                    vmin=-vlims[TF_type][chan],
                    vmax=vlims[TF_type][chan],
                    cmap='seismic'
                )

                ax.vlines(
                    x=int(stretch_point_TF / 2),
                    ymin=frex[0],
                    ymax=frex[-1],
                    color="k",
                    linewidth=2
                )

                ax.set_yscale('log')

                ax.set_yticks(ticks_freq)
                ax.set_yticklabels([str(t) for t in ticks_freq])

                if row_i == 0:
                    ax.set_title(
                        f"{cond} {xr_dict_data[cond].shape[1]}"
                    )

                if col_i == 0:
                    ax.set_ylabel(TF_type)

            # one colorbar for the entire row
            cbar = fig.colorbar(
                pcm,
                ax=axs[row_i, :],
                orientation="vertical",
                fraction=0.02,
                pad=0.02
            )

        plt.suptitle(f"{sujet} {chan}")

        # plt.show()

        fig.savefig(os.path.join(path_export_TF, f"{sujet}_{chan}_TF.jpeg"))






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        tf_relative_patientwise(sujet)

    

        


    























