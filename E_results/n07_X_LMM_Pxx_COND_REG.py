



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







########################################
######## GENERATE ALLDF COND ########
########################################

def get_LMM_COND_df_Pxx_patient_wise_VCCP():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['A+', 'R-', 'C-']

    path_import_res_LMM = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx', 'LMM_res', 'COND')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_LMM, f"RES_VCCP_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx





def get_LMM_COND_df_Pxx_patient_wise_CB():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['condCB']

    path_import_res_LMM = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx', 'LMM_res', 'COND')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_LMM, f"RES_attention_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx





def get_LMM_COND_df_Pxx_patient_wise_HV():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['condHV']

    path_import_res_LMM = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx', 'LMM_res', 'COND')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_LMM, f"RES_HV_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx







########################################
######## GENERATE ALLDF REG ########
########################################




def get_LMM_REG_df_Pxx_patient_wise_VCCP():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['A_scaled', 'R_scaled', 'C_scaled', 'A_scaled:R_scaled', 'A_scaled:C_scaled', 'R_scaled:C_scaled', 'A_scaled:R_scaled:C_scaled']

    path_import_res_LMM = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx', 'LMM_res', 'REG')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_LMM, f"RES_VCCP_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx






def get_LMM_REG_df_Pxx_patient_wise_CB():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['A_scaled', 'R_scaled', 'condCB', 'A_scaled:R_scaled', 'A_scaled:condCB', 'R_scaled:condCB', 'A_scaled:R_scaled:condCB']

    path_import_res_LMM = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx', 'LMM_res', 'REG')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_LMM, f"RES_CBCOND_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx





def get_LMM_REG_df_Pxx_patient_wise_CO2():

    phase_cycle_list = ['inspi', 'expi']
    factor_select = ['A_scaled', 'R_scaled', 'C-', 'A_scaled:R_scaled', 'A_scaled:C-', 'R_scaled:C-', 'A_scaled:R_scaled:C-']

    path_import_res_LMM = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx', 'LMM_res', 'REG')

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(os.path.join(path_import_res_LMM, f"RES_VCCP_CO2COND_{band}_{phase_cycle}_{ROI}_Pxx.xlsx"))
                _df = _df.query(f"term in {factor_select}").reset_index(drop=True)
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'}).reset_index(drop=True)

    return df_R_Pxx

















########################
######## PXX ########
########################




def export_res_Pxx_patient_wise(sujet):

    #### config
    phase_cycle_list = ['inspi', 'expi']
    path_export_res = os.path.join(path_results, 'Pxx', 'LMM')

    #### load df
    df_R_Pxx_COND_VCCP = get_LMM_COND_df_Pxx_patient_wise_VCCP()
    df_R_Pxx_COND_CB = get_LMM_COND_df_Pxx_patient_wise_CB()
    df_R_Pxx_COND_HV = get_LMM_COND_df_Pxx_patient_wise_HV()

    df_R_Pxx_VCCP_REG = get_LMM_REG_df_Pxx_patient_wise_VCCP()
    df_R_Pxx_CB_REG = get_LMM_REG_df_Pxx_patient_wise_CB()
    df_R_Pxx_CO2_REG = get_LMM_REG_df_Pxx_patient_wise_CO2()

    chanlist, chanlist_aux, chanlist_all, localist = get_chanlist_localist(sujet)

    label_name_nplot = {}
    
    for ROI in ROI_short_list:

        _nplot = localist.query(f"loca == '{ROI}'").size
        label_name_nplot[ROI] = f"{ROI} c({_nplot})"

    for ROI in ROI_short_list:
        df_R_Pxx_COND_VCCP['ROI'] = df_R_Pxx_COND_VCCP['ROI'].replace({ROI: label_name_nplot[ROI]})
        df_R_Pxx_COND_CB['ROI'] = df_R_Pxx_COND_CB['ROI'].replace({ROI: label_name_nplot[ROI]})
        df_R_Pxx_COND_HV['ROI'] = df_R_Pxx_COND_HV['ROI'].replace({ROI: label_name_nplot[ROI]})

        df_R_Pxx_VCCP_REG['ROI'] = df_R_Pxx_VCCP_REG['ROI'].replace({ROI: label_name_nplot[ROI]})
        df_R_Pxx_CB_REG['ROI'] = df_R_Pxx_CB_REG['ROI'].replace({ROI: label_name_nplot[ROI]})
        df_R_Pxx_CO2_REG['ROI'] = df_R_Pxx_CO2_REG['ROI'].replace({ROI: label_name_nplot[ROI]})

    ROI_plot_short_list = df_R_Pxx_COND_VCCP['ROI'].unique().tolist()

    #### config figure
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    #### COND
    protocol_dict = {'VCCP_COND' : df_R_Pxx_COND_VCCP, 'CB_COND' : df_R_Pxx_COND_CB, 'HV_COND' : df_R_Pxx_COND_HV}
    LMM_param_dict = {'VCCP_COND' : ['A+', 'R-', 'C-'], 'CB_COND' : ['condCB'], 'HV_COND' : ['condHV']}

    for protocol, df_protocol in protocol_dict.items():
    
        #band = 'theta'
        for band in freq_band_dict:

            fig = make_subplots(
                rows=len(LMM_param_dict[protocol]),
                cols=1,
                shared_xaxes=True,
                subplot_titles=LMM_param_dict[protocol]
            )

            y_max = df_protocol.query(f"band == '{band}' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

            for r, LMM_param in enumerate(LMM_param_dict[protocol], start=1):

                _df_plot = df_protocol.query(f"band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}").copy()

                _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                for phase_cycle in _df_plot["phase_cycle"].unique():

                    df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                    fig.add_bar(
                        x=df_term["ROI"],
                        y=df_term["estimate"],
                        name=phase_cycle,
                        text=df_term["sig"],
                        textposition="outside",
                        marker_color=color_map[phase_cycle],
                        row=r,
                        col=1
                    )

                    margin = y_max * 0.50   # 15% extra space

                    fig.update_yaxes(
                        range=[-y_max - margin, y_max + margin],
                        row=r,
                        col=1
                    )

                    # fig.update_xaxes(
                    #     categoryorder="array",
                    #     categoryarray=ROI_order,
                    #     row=r,
                    #     col=1
                    # )

            fig.update_layout(
                title=f"COND {protocol} {band} all ROI",
                template="simple_white",
                barmode="group",
                xaxis_tickangle=-45,
                yaxis_title="estimate",
                legend_title="phase_cycle",
                height=320 * len(phase_cycle_list),  # scale height
                width=500,
            )

            # fig.show()

            fig.write_html(
                os.path.join(path_export_res, protocol, 'patientwise', f"{protocol}_{band}_threshROI_LMM.html"),
                include_plotlyjs="cdn"
            )



    #### REG    
    protocol_dict = {'VCCP_REG' : df_R_Pxx_VCCP_REG, 'CB_REG' : df_R_Pxx_CB_REG, 'CO2_REG' : df_R_Pxx_CO2_REG}
    LMM_param_dict = {'VCCP_REG' : ['A_scaled', 'R_scaled', 'C_scaled', 'A_scaled:R_scaled', 'A_scaled:C_scaled', 'R_scaled:C_scaled', 'A_scaled:R_scaled:C_scaled'], 
                      'CB_REG' : ['A_scaled', 'R_scaled', 'condCB', 'A_scaled:R_scaled', 'A_scaled:condCB', 'R_scaled:condCB', 'A_scaled:R_scaled:condCB'], 
                      'CO2_REG' : ['A_scaled', 'R_scaled', 'C-', 'A_scaled:R_scaled', 'A_scaled:C-', 'R_scaled:C-', 'A_scaled:R_scaled:C-']}

    for protocol, df_protocol in protocol_dict.items():
    
        #band = 'theta'
        for band in freq_band_dict:

            fig = make_subplots(
                rows=len(LMM_param_dict[protocol]),
                cols=1,
                shared_xaxes=True,
                subplot_titles=LMM_param_dict[protocol]
            )

            y_max = df_protocol.query(f"band == '{band}' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

            for r, LMM_param in enumerate(LMM_param_dict[protocol], start=1):

                _df_plot = df_protocol.query(f"band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}").copy()

                _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                for phase_cycle in _df_plot["phase_cycle"].unique():

                    df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                    fig.add_bar(
                        x=df_term["ROI"],
                        y=df_term["estimate"],
                        name=phase_cycle,
                        text=df_term["sig"],
                        textposition="outside",
                        marker_color=color_map[phase_cycle],
                        row=r,
                        col=1
                    )

                    # margin = y_max * 0.50   # 15% extra space

                    # fig.update_yaxes(
                    #     range=[-y_max - margin, y_max + margin],
                    #     row=r,
                    #     col=1
                    # )

                    # fig.update_xaxes(
                    #     categoryorder="array",
                    #     categoryarray=ROI_order,
                    #     row=r,
                    #     col=1
                    # )

            fig.update_layout(
                title=f"REG {protocol} {band} all ROI",
                template="simple_white",
                barmode="group",
                xaxis_tickangle=-45,
                yaxis_title="estimate",
                legend_title="phase_cycle",
                height=200 * len(LMM_param_dict[protocol]),  # scale height
                width=500,
            )

            # fig.show()

            fig.write_html(
                os.path.join(path_export_res, protocol, 'patientwise', f"{protocol}_{band}_threshROI_LMM.html"),
                include_plotlyjs="cdn"
            )











################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        export_res_Pxx_patient_wise(sujet)

    

        


    























