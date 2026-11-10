



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *






################################
######## TF RELATIVE ########
################################



def tf_relative_patientwise(sujet, ref_type):

    print(f"{sujet} {ref_type}")

    #### config
    path_load_TF = os.path.join(path_precompute, 'TF')
    path_export_TF = os.path.join(path_results, 'TF', sujet)
        
    TF_type_list = ['absolute', 'diff']

    #### load extract vlim
    respfeatures = get_respfeatures(sujet)
    chanlist, localist = get_chanlocalist(sujet, ref_type)

    xr_dict_data = {cond : xr.open_dataarray(os.path.join(path_load_TF, f"{sujet}_{cond}_tf_stretch_{ref_type}.nc")) for cond in conditions}
    
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

        fig.savefig(os.path.join(path_export_TF, f"{sujet}_{chan}_TF_{ref_type}.jpeg"))

        plt.close('all')





def tf_relative_population(ref_type):

    #### config
    path_export_fig = os.path.join(path_results, 'TF', 'population')
    path_load = os.path.join(path_precompute, 'TF')

    TF_type_list = ['absolute', 'diff']

    df_loca_allsujet = get_df_loca_allsujet(ref_type)

    ROI_list = df_loca_allsujet['loca_modified'].unique()

    #### ROI computation
    for ROI in ROI_list:

        print(ROI)

        sujet_list_ROI = df_loca_allsujet.query(f"loca_modified == '{ROI}'")['sujet'].unique()

        #### get sujet extraction params
        ROI_TF = {}

        for cond in conditions:

            ROI_TF[cond] = []
        
            for sujet in sujet_list_ROI:

                xr_tf_data = xr.open_dataarray(os.path.join(path_load, f"{sujet}_{cond}_tf_stretch_{ref_type}.nc"))
                chanlist_extraction = df_loca_allsujet.query(f"loca_modified == '{ROI}' and sujet == '{sujet}'")['chan'].values
                xr_tf_data_chansel = xr_tf_data.sel(chan=chanlist_extraction)
                tf_median = xr_tf_data_chansel.median('cycle').median('chan').values

                ROI_TF[cond].append(tf_median.reshape(1,tf_median.shape[0],tf_median.shape[1]))

            ROI_TF[cond] = np.concat(ROI_TF[cond], axis=0)

        #### extract vlims
        min = np.array([ROI_erp[cond].min() for cond in conditions]).min()
        max = np.array([ROI_erp[cond].max() for cond in conditions]).max()

        #### get vlims
        vlims = {}
        
        for TF_type in TF_type_list:
    
            vlims[TF_type] = {}
        
            vlim = []

            data_baseline = np.median(ROI_TF['RB'], axis=0)
            vlim.append([data_baseline.min(), data_baseline.max()])

            for cond in conditions_oc:

                if TF_type == 'absolute':
                    
                    vlim.append([np.median(ROI_TF[cond], axis=0).min(), np.median(ROI_TF[cond], axis=0).max()])

                elif TF_type == 'relative':

                    _tf_diff = np.median((ROI_TF[cond] - ROI_TF['RB']), axis=0)
                    vlim.append([_tf_diff.min(), _tf_diff.max()])

            vlims[TF_type] = np.abs(np.array(vlim)).max()

        #### plot
        phase_vec = np.arange(stretch_point_TF)
        ticks_freq = [2, 8, 12, 30, 50, 60, 150]

        fig, axs = plt.subplots(
            ncols=len(conditions),
            nrows=len(TF_type_list),
            figsize=(16, 10)
        )

        for row_i, TF_type in enumerate(TF_type_list):

            for col_i, cond in enumerate(conditions):

                if cond == 'RB':
                    tf_plot = np.median(ROI_TF['RB'], axis=0)

                else:
                    if TF_type == 'absolute':
                        tf_plot = np.median(ROI_TF[cond], axis=0)

                    elif TF_type == 'diff':
                        tf_plot = np.median((ROI_TF[cond] - ROI_TF['RB']), axis=0)

                ax = axs[row_i, col_i]

                pcm = ax.pcolormesh(
                    phase_vec,
                    frex,
                    tf_plot,
                    vmin=-vlims[TF_type],
                    vmax=vlims[TF_type],
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
                        f"{cond}"
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

        plt.suptitle(f"ROI:{ROI} {ROI_TF[cond].shape[0]}")

        # plt.show()

        fig.savefig(os.path.join(path_export_TF, 'population', f"{ROI}_TF_{ref_type}.jpeg"))

        plt.close('all')







################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #ref_type = 'mn'
    for ref_type in ref_type_list:

        #sujet = sujet_list[0]
        for sujet in sujet_list:

            tf_relative_patientwise(sujet, ref_type)

    

        


    























