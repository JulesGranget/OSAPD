



from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *


import plotly.express as px
from plotly.subplots import make_subplots






################################
######## STATS FUNCTIONS ########
################################

def p_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''





########################
######## PXX ########
########################




def export_res_Pxx_general():

    #### config
    sujet_thresh = 3
    phase_cycle_list = ['inspi', 'expi']

    #### load df
    df_loca_allsujet = get_df_loca_allsujet(ref_type)
    ROI_list = df_loca_allsujet['loca_modified'].unique()

    phase_cycle_list = ['inspi', 'expi']
    
    path_load_LMM_res = os.path.join(path_precompute, 'Pxx', 'LMM', 'df_stats')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_list:

                _df = pd.read_excel(os.path.join(path_load_LMM_res, f"PXX_COND_OC_RES_{band}_{phase_cycle}_{ROI}_LMM.xlsx"))
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value', 'std.error']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue', 'std.error' : 'std_error'})
    
    #### anatomy
    label_name_nsujet = {}
    
    for ROI in ROI_list:

        _nsujet = df_loca_allsujet.query(f"loca_modified == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} (p{_nsujet})"

    for ROI in ROI_list:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca_modified').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_list:
        df_ROI_scount['loca_modified'] = df_ROI_scount['loca_modified'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca_modified'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    #### main plot
        #### inspi/expi
    df_main_plot = df_R_Pxx.copy()
    df_main_plot['term'] = df_main_plot['term'].replace({'oc_ratio' : 'OR', 'condoc300' : 'O300', 'condoc600' : 'O600', 
                                                         'oc_ratio:condoc300' : 'OR:O300', 'oc_ratio:condoc600' : 'OR:O600'})

    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['OR', 'O300', 'O600', 'OR:O300', 'OR:O600']
    
    ROI_order = ['ant (p12)', 'post (p11)', 'STN (p7)']

    #band = 'theta'
    for band in freq_band_dict:

        fig = make_subplots(
            rows=1,
            cols=len(LMM_param_list),
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        for c, LMM_param in enumerate(LMM_param_list, start=1):

            _df_plot = df_main_plot.query(
                f"band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}"
            ).copy()

            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            phases = list(_df_plot["phase_cycle"].unique())

            # keep phase/hue points close together within each ROI
            phase_offsets = dict(zip(phases, np.linspace(-0.12, 0.12, len(phases))))

            # ROI positions remain clearly separated
            roi_to_x = {roi: i for i, roi in enumerate(ROI_order)}

            y_min = (_df_plot["estimate"] - _df_plot["std_error"]).min()
            y_max = (_df_plot["estimate"] + _df_plot["std_error"]).max()

            y_range = y_max - y_min
            star_offset = 0.04 * y_range
            margin = 0.15 * y_range

            for phase_cycle in phases:

                df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle].copy()

                x_pos = np.array([roi_to_x[roi] for roi in df_term["ROI"]]) + phase_offsets[phase_cycle]

                # estimate + SE
                fig.add_scatter(
                    x=x_pos,
                    y=df_term["estimate"],
                    mode="markers",
                    name=phase_cycle,
                    marker=dict(
                        color=color_map[phase_cycle],
                        size=10
                    ),
                    error_y=dict(
                        type="data",
                        array=df_term["std_error"],
                        visible=True,
                        color=color_map[phase_cycle],
                        thickness=1.5,
                        width=4
                    ),
                    showlegend=(c == 1),
                    row=1,
                    col=c
                )

                # significance stars above error bar
                fig.add_scatter(
                    x=x_pos,
                    y=df_term["estimate"] + df_term["std_error"] + star_offset,
                    mode="text",
                    text=df_term["sig"],
                    textposition="top center",
                    textfont=dict(size=16),
                    showlegend=False,
                    hoverinfo="skip",
                    row=1,
                    col=c
                )

            fig.update_yaxes(
                range=[y_min - margin, y_max + 2 * margin],
                zeroline=True,
                zerolinecolor="grey",
                zerolinewidth=1,
                row=1,
                col=c
            )

            fig.update_xaxes(
                tickmode="array",
                tickvals=list(range(len(ROI_order))),
                ticktext=ROI_order,
                tickangle=-45,
                row=1,
                col=c
            )

        fig.update_layout(
            title=f"{band}",
            template="simple_white",
            yaxis_title="<i>β</i> estimate",
            legend_title="phase_cycle",
            height=500,
            width=300 * len(LMM_param_list)
        )

        # fig.show()

        outdir = os.path.join(path_results, "Pxx", "LMM")
        fig.write_html(
            os.path.join(outdir, f"{band}_Pxx_COND_OC_LMM.html"),
            include_plotlyjs="cdn"
        )

        if band in ['theta', 'beta', 'gamma']:

            fig.write_image(
                os.path.join(path_paper_figure_export, f"fig00_{band}_threshROI_LMM.svg")
            )

    #### example plot
    df_example_plot = df_R_Pxx.copy()
    df_example_plot['term'] = df_example_plot['term'].replace({'oc_ratio' : 'OR', 'condoc300' : 'O300', 'condoc600' : 'O600', 
                                                             'oc_ratio:condoc300' : 'OR:O300', 'oc_ratio:condoc600' : 'OR:O600'})

    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['OR', 'O300', 'O600', 'OR:O300', 'OR:O600']
    
    ROI_sel = 'STN (p7)'
    band_sel = 'theta'

    fig = make_subplots(
        rows=1,
        cols=len(LMM_param_list),
        shared_xaxes=False,
        subplot_titles=LMM_param_list,
        horizontal_spacing=0.1
    )

    _df_max = df_example_plot.query(
        f"band == '{band_sel}' and term != '(Intercept)' and phase_cycle != 'whole' and ROI == '{ROI_sel}'"
    ).copy()

    y_min = (_df_max["estimate"] - _df_max["std_error"]).min()
    y_max = (_df_max["estimate"] + _df_max["std_error"]).max()
    y_range = y_max - y_min if y_max != y_min else 1

    star_offset = 0.04 * y_range
    margin = 0.15 * y_range

    for c, LMM_param in enumerate(LMM_param_list, start=1):

        _df_plot = df_example_plot.query(
            f"band == '{band_sel}' and term == '{LMM_param}' and phase_cycle != 'whole' and ROI == '{ROI_sel}'"
        ).copy()

        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        phases = list(_df_plot["phase_cycle"].unique())
        _offset = 0.15
        phase_offsets = {'inspi' : -_offset, 'expi' : _offset}

        for phase_cycle in phases:

            df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle].copy()
            x_pos = np.zeros(len(df_term)) + phase_offsets[phase_cycle]

            fig.add_scatter(
                x=x_pos,
                y=df_term["estimate"],
                mode="markers",
                name=phase_cycle,
                marker=dict(color=color_map[phase_cycle], size=10),
                error_y=dict(
                    type="data",
                    array=df_term["std_error"],
                    visible=True,
                    color=color_map[phase_cycle],
                    thickness=1.5,
                    width=4
                ),
                showlegend=(c == 1),
                row=1,
                col=c
            )

            fig.add_scatter(
                x=x_pos,
                y=df_term["estimate"] + df_term["std_error"] + star_offset,
                mode="text",
                text=df_term["sig"],
                textposition="top center",
                textfont=dict(size=16),
                showlegend=False,
                hoverinfo="skip",
                row=1,
                col=c
            )

        fig.update_yaxes(
            range=[y_min - margin, y_max + 2 * margin],
            tickfont=dict(size=14),
            title_text="β estimate" if c == 1 else None,
            zeroline=True,
            zerolinecolor="grey",
            zerolinewidth=1,
            row=1,
            col=c
        )

        fig.update_xaxes(
            range=[-1, 1],
            tickmode="array",
            tickvals=[0],
            ticktext=[ROI_sel],
            tickfont=dict(size=14),
            row=1,
            col=c
        )

    fig.update_layout(
        title=f"{band_sel} {ROI_sel}",
        template="simple_white",
        legend_title="phase_cycle",
        height=500,
        width=200 * len(LMM_param_list)
    )

    # fig.show()

    fig.write_image(
        os.path.join(path_paper_figure_export, f"fig00_EXAMPLE_{band_sel}_Pxx_LMM.svg")
    )

    #### export df data
    df_export = df_R_Pxx.query(f"phase_cycle != 'whole' and ROI in {ROI_plot_short_list} and term != '(Intercept)'")
    
    output_df = os.path.join(path_results, "Pxx", "LMM")
    filename = os.path.join(output_df, "df_Pxx_LMM.xlsx")
    df_export.to_excel(filename)












################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_res_Pxx_general()
    

        


    























