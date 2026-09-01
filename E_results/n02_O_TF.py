



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *






################################
######## TF RELATIVE ########
################################



def tf_relative_patientwise(sujet):

    #### config
    path_load_TF = os.path.join(path_precompute, 'TF', sujet)
    path_export_TF = os.path.join(path_results, 'TF', 'patientwise', sujet)

    #### load    
    xr_allcond_sujet = xr.open_dataarray(os.path.join(path_load_TF, f"{sujet}_ALLCOND_tf_stretch_RELABEL.nc"))
    xr_allcond_sujet = xr_allcond_sujet.where(xr_allcond_sujet.loca.isin(ROI_short_list), drop=True)
    
    xr_VCCP = xr_allcond_sujet.where(~xr_allcond_sujet.cond.isin(['RB', 'HV']), drop=True)
    xr_VCCP_rel = xr_VCCP.where(~xr_VCCP.cond.isin(['AoRoCo']), drop=True) - xr_VCCP.sel(cond='AoRoCo') 
    
    xr_HV = xr_allcond_sujet.where(xr_allcond_sujet.cond.isin(['RB', 'HV']), drop=True)

    xr_CB = xr_allcond_sujet.where(xr_allcond_sujet.cond.isin(['RB', 'AoRoCo']), drop=True)

    #### plot
    protocol_dict = {'VCCP' : xr_VCCP_rel, 'HV' : xr_HV, 'CB' : xr_CB}
    cond_plot_protocol = {'VCCP' : [cond for cond in conditions if cond not in ['RB', 'HV', 'AoRoCo']], 
                          'HV' : ['RB', 'HV'],
                          'CB' : ['RB', 'AoRoCo']}
    fig_size_list = {'VCCP' : (24,5), 'HV' : (10,5), 'CB' : (10,5)}
    phase_vec = xr_allcond_sujet.phase.data
    ticks_freq = [2, 8, 12, 30, 50, 60, 150]

    for protocol, xr_protocol in protocol_dict.items():

        for ROI in ROI_short_list:

            #### sel
            nchan = xr_protocol.where(xr_protocol.loca == ROI, drop=True).shape[1]
            xr_plot = xr_protocol.where(xr_protocol.loca == ROI, drop=True).median('chan')
            cond_sel_plot = cond_plot_protocol[protocol]
            
            #### vlim
            _vlim = np.abs([np.percentile(xr_plot, percentile_plot[0]), np.percentile(xr_plot, percentile_plot[1])]).max()

            fig, axs = plt.subplots(ncols=len(cond_sel_plot), figsize=fig_size_list[protocol])
            
            for cond_i, cond in enumerate(cond_sel_plot):
            
                ax = axs[cond_i]
                pcm = ax.pcolormesh(phase_vec, frex, xr_plot.sel(cond=cond), vmin=-_vlim, vmax=_vlim, cmap='seismic')
                ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
                ax.set_yscale('log')
                
                ax.set_yticks(ticks_freq)
                ax.set_yticklabels([str(t) for t in ticks_freq])

                ax.set_title(f"{cond}")

            cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
            plt.suptitle(f"{protocol} {sujet} {ROI} c({nchan})")

            # plt.show()

            fig.savefig(os.path.join(path_export_TF, f"{sujet}_{protocol}_{ROI}.jpg"))

        plt.close('all')






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        tf_relative_patientwise(sujet)

    

        


    























