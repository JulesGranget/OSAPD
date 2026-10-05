


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import matplotlib.colors as mcolors


################################
######## ERP PLOT ########
################################



def plot_ERP_res(sujet, ref_type):

    #### config
    path_export_fig = os.path.join(path_results, 'ERP', sujet)

    freq_cond_list = ['WB', 'HFA']
    conditions_sel = ['RB', 'oc300', 'oc600']
    time_vec_ERP = np.arange(-chunk_time_prepost, chunk_time_prepost, 1/srate)
    time_vec_plot = (-1 < time_vec_ERP) & (time_vec_ERP < 2)
    path_load = os.path.join(path_precompute, 'ERP', sujet)
    chanlist, localist = get_chanlocalist(sujet, ref_type)

    ### load ERP
    data_ERP = {}

    for freq_type in freq_cond_list:

        data_ERP[freq_type] = {}

        for cond in conditions_sel:        

            data_ERP[freq_type][cond] = xr.open_dataarray(os.path.join(path_load, f"{sujet}_{cond}_{freq_type}_erp_{ref_type}.nc"))

    ### load stats
    stats_ERP = {}

    for freq_type in freq_cond_list:

        stats_ERP[freq_type] = xr.open_dataarray(os.path.join(path_load, f"{sujet}_perm_stats_{freq_type}_{ref_type}.nc"))

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
                # _min.append(_sig_mSE.min())
                # _max.append(_sig_pSE.max())
                _min.append(_sig.min()*1.1)
                _max.append(_sig.max()*1.1)

                SE[freq_type][chan][cond] = [_sig_mSE, _sig_pSE]

            vlims[freq_type][chan] = [np.array(_min).min(), np.array(_max).max()]

    if debug:

        cond_color = {"RB": "blue", "oc300": "orange", "oc600": "red"}

        chan = chanlist[0]

        fig, ax = plt.subplots()
        
        for cond in conditions_sel:

            _sig_mat = data_ERP[freq_type][cond].sel(chan=chan).data
            _sig_med = np.median(_sig_mat, axis=0)
            _sig_mSE, _sig_pSE = SE[freq_type][chan][cond][0], SE[freq_type][chan][cond][1]

            ax.plot(time_vec_ERP[time_vec_plot], _sig_med[time_vec_plot], color=cond_color[cond], label=f"{cond}:c{_sig_mat.shape[0]}")
            ax.fill_between(time_vec_ERP[time_vec_plot], _sig_mSE[time_vec_plot], _sig_pSE[time_vec_plot], color=cond_color[cond], alpha=0.5)

        plt.legend()
        fig.show()

        fig, ax = plt.subplots()
                
        cond = conditions_sel[0]

        _sig_mat = data_ERP[freq_type][cond].sel(chan=chan).data
        _sig_med = np.median(_sig_mat, axis=0)

        for cycle_i in range(_sig_mat.shape[0]):
            ax.plot(time_vec_ERP[time_vec_plot], _sig_mat[cycle_i,time_vec_plot], alpha=0.2, label=f"{cond}:c{_sig_mat.shape[0]}")
        ax.plot(time_vec_ERP[time_vec_plot], _sig_med[time_vec_plot], color=cond_color[cond], label=f"{cond}:c{_sig_mat.shape[0]}")

        fig.show()

                
    #### plot
    # cond_color = {'RB' : 'g', 'oc300' : 'k', 'oc600' : 'r'}
    cond_color = {"RB": "blue", "oc300": "orange", "oc600": "red"}

    #chan_i, chan = 0, chanlist[0]
    for chan_i, chan in enumerate(chanlist):

        #### figure
        def color_to_rgba(color, alpha=0.15):
            r, g, b, _ = mcolors.to_rgba(color)
            return f"rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, {alpha})"

        fig = make_subplots(
            rows=len(freq_cond_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=[str(freq_type) for freq_type in freq_cond_list],
            vertical_spacing=0.06
        )

        for row_i, freq_type in enumerate(freq_cond_list, start=1):

            for cond in conditions_sel:

                _sig_mat = data_ERP[freq_type][cond].sel(chan=chan).data
                _sig_med = np.median(_sig_mat, axis=0)
                _sig_mSE, _sig_pSE = SE[freq_type][chan][cond][0], SE[freq_type][chan][cond][1]

                x = time_vec_ERP[time_vec_plot]
                y_med = _sig_med[time_vec_plot]

                ymin, ymax = vlims[freq_type][chan]

                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y_med,
                        mode="lines",
                        line=dict(color=cond_color[cond], width=2),
                        name=f"{cond}:c{_sig_mat.shape[0]}",
                        legendgroup=cond,
                        showlegend=(row_i == 1)
                    ),
                    row=row_i,
                    col=1
                )

                #### significance
                if cond != 'RB':
                    sig_bool = stats_ERP[freq_type].sel(cond=cond, chan=chan).data.astype('bool')
                    sig_bool = sig_bool[time_vec_plot]
                    # sig_bool[int(sig_bool.size/2):] = 1
                    x_sig = x[sig_bool]

                    y_sig_val = ymin + 0.05 * (ymax - ymin)
                    y_sig = np.full(x_sig.size, y_sig_val)

                    fig.add_trace(
                        go.Scatter(
                            x=x_sig,
                            y=y_sig,
                            mode="markers",
                            marker=dict(size=5, color="black"),
                            showlegend=False,
                            hoverinfo="skip"
                        ),
                        row=row_i,
                        col=1
                    )

            fig.update_yaxes(range=[ymax, ymin], row=row_i, col=1)
            fig.add_vline(x=0, line_width=1.5, line_color="black", row=row_i, col=1)

        fig.update_xaxes(title_text="Time (s)", row=len(freq_cond_list), col=1)
        fig.update_xaxes(tickfont=dict(size=12))
        fig.update_yaxes(tickfont=dict(size=12))

        fig.update_layout(
            title=dict(text=f"{sujet} {chan} {ref_type}", x=0.5, xanchor="center"),
            template="simple_white",
            height=300 * len(freq_cond_list),
            width=1000,
            hovermode="x unified",
            legend=dict(title="Condition", orientation="v"),
            margin=dict(l=80, r=80, t=100, b=70)
        )

        # fig.show()

        fig.write_html(os.path.join(path_export_fig, f"{sujet}_{chan}_{ref_type}.html"))

        ####
        # fig, axs = plt.subplots(nrows=len(freq_cond_list), ncols=1)
        
        # plt.suptitle(f'{chan}')

        # fig.set_figheight(10)
        # fig.set_figwidth(10)

        # #row_i, freq_type = 0, freq_cond_list[0]
        # for row_i, freq_type in enumerate(freq_cond_list):

        #     ax = axs[row_i]

        #     #cond = conditions_sel[0]
        #     for cond in conditions_sel:    

        #         _sig_mat = data_ERP[freq_type][cond].sel(chan=chan).data
        #         _sig_med = np.median(_sig_mat, axis=0)
        #         _sig_mSE, _sig_pSE = SE[freq_type][chan][cond][0], SE[freq_type][chan][cond][1]

        #         ax.plot(time_vec_ERP[time_vec_plot], _sig_med[time_vec_plot], color=cond_color[cond], label=f"{cond}:c{_sig_mat.shape[0]}")
        #         ax.fill_between(time_vec_ERP[time_vec_plot], _sig_mSE[time_vec_plot], _sig_pSE[time_vec_plot], color=cond_color[cond], alpha=0.5)

        #     ax.set_title(f"{freq_type}", fontweight='bold')

        #     ax.set_ylim(vlims[freq_type][chan][0], vlims[freq_type][chan][1])

        #     ax.invert_yaxis()

        #     ax.vlines(0, ymin=vlims[freq_type][chan][0], ymax=vlims[freq_type][chan][1], colors='k')

        #     ax.legend()

        # fig.tight_layout()

        # # plt.show()

        # #### save
        # fig.savefig(os.path.join(path_export_fig, f'{sujet}_{chan}.jpeg'), dpi=150)

        # fig.clf()
        # plt.close('all')


            






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        #ref_type = 'mn'
        for ref_type in ref_type_list:

            plot_ERP_res(sujet, ref_type)

    
