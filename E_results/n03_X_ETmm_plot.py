



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







################################
######## MEDIAN SIG ########
################################



def plot_median_ETmm_allcond(sujet):

    #### config
    path_export_res = os.path.join(path_results, 'eyet', 'mean')
    path_load_sig = os.path.join(path_precompute, 'ET')

    #### load df reg
    xr_ETmm = xr.load_dataarray(os.path.join(path_load_sig, f"{sujet}_ET_allcond.nc"))

    #### sel xr
    xr_ETmm_VCCP = xr_ETmm.where(~xr_ETmm.cond.isin(['RB', 'HV']), drop=True)
    xr_ETmm_HV = xr_ETmm.where(xr_ETmm.cond.isin(['RB', 'HV']), drop=True)
    xr_ETmm_CB = xr_ETmm.where(xr_ETmm.cond.isin(['RB', 'AoRoCo']), drop=True)
    
    #### VCCP
    xr_plot = xr_ETmm_VCCP.median('cycle').where(~xr_ETmm_VCCP.median('cycle').cond.isin(['AoRoCo']), drop=True) - xr_ETmm_VCCP.median('cycle').sel(cond='AoRoCo')
    
    fig, ax = plt.subplots(figsize=(8, 5))

    xr_plot.plot.line(x='phase', hue='cond', ax=ax)

    ax.set_title('ETmm relative AoRoCo')

    # plt.show()

    fig.savefig(os.path.join(path_export_res, f"{sujet}_VCCP_ALLCOND_median_ETmm.png"))

    #### HV
    xr_plot = xr_ETmm_HV.median('cycle')
    
    fig, ax = plt.subplots(figsize=(8, 5))

    xr_plot.plot.line(x='phase', hue='cond', ax=ax)

    ax.set_title('ETmm HV')

    # plt.show()

    fig.savefig(os.path.join(path_export_res, f"{sujet}_HV_ALLCOND_median_ETmm.png"))

    #### CB
    xr_plot = xr_ETmm_CB.median('cycle')
    
    fig, ax = plt.subplots(figsize=(8, 5))

    xr_plot.plot.line(x='phase', hue='cond', ax=ax)

    ax.set_title('ETmm CB')

    # plt.show()

    fig.savefig(os.path.join(path_export_res, f"{sujet}_CB_ALLCOND_median_ETmm.png"))








########################
######## REG ########
########################




def plot_REG_patient_wise(sujet):

    #### config
    path_export_res = os.path.join(path_results, 'eyet', 'REG', 'patientwise')

    #### load df reg
    path_load_reg = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'ETmm')

    df_reg_VCCP = pd.read_excel(os.path.join(path_load_reg, f"{sujet}_df_reg_ETmm_VCCP_R.xlsx")).drop(columns=['Unnamed: 0'])
    df_reg_HV = pd.read_excel(os.path.join(path_load_reg, f"{sujet}_df_reg_ETmm_HV_R.xlsx")).drop(columns=['Unnamed: 0'])
    df_reg_CB = pd.read_excel(os.path.join(path_load_reg, f"{sujet}_df_reg_ETmm_CB_R.xlsx")).drop(columns=['Unnamed: 0'])

    #### plot
    cond_list_reg = {'VCCP' : df_reg_VCCP, 'HV' : df_reg_HV, 'CB' : df_reg_CB}

    for _cond_protocol, df_protocol in cond_list_reg.items():

        df_plot = df_protocol
        col_melt = ['A_scaled', 'R_scaled', 'C_scaled']
        df_plot = pd.melt(df_plot, id_vars=[col for col in df_plot.columns if col not in col_melt], value_vars=col_melt, var_name='factor', value_name='rscore')

        fig = sns.lmplot(data=df_plot, x='rscore', y='ETmm_scaled', col='factor')
        fig.fig.subplots_adjust(top=0.92, wspace=0.35, hspace=0.45)
        fig.fig.suptitle(f"{_cond_protocol}")
        
        # fig.fig.show()

        fig.fig.savefig(os.path.join(path_export_res, f"{sujet}_{_cond_protocol}_REG.png"))




def plot_STATS_REG_patient_wise(sujet):

    #### config
    path_export_res = os.path.join(path_results, 'eyet', 'REG', 'patientwise')

    #### load df STATS reg
    path_load_reg_stats = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'ETmm', 'LMM_res', 'REG')

    df_stats_VCCP = pd.read_excel(os.path.join(path_load_reg_stats, f"RES_{sujet}_VCCP_ETmm.xlsx"))
    df_stats_HV = pd.read_excel(os.path.join(path_load_reg_stats, f"RES_{sujet}_HV_ETmm.xlsx"))
    df_stats_CB = pd.read_excel(os.path.join(path_load_reg_stats, f"RES_{sujet}_CB_ETmm.xlsx"))

    #### REG    
    protocol_dict = {'VCCP_REG' : df_stats_VCCP, 'HV_REG' : df_stats_HV, 'CB_REG' : df_stats_CB}
    LMM_param_dict = {'VCCP_REG' : ['A_scaled', 'R_scaled', 'C_scaled', 'A_scaled:R_scaled', 'A_scaled:C_scaled', 'R_scaled:C_scaled', 'A_scaled:R_scaled:C_scaled'], 
                      'HV_REG' : ['A_scaled', 'R_scaled', 'C_scaled', 'A_scaled:R_scaled', 'A_scaled:C_scaled', 'R_scaled:C_scaled', 'A_scaled:R_scaled:C_scaled'], 
                      'CB_REG' : ['A_scaled', 'R_scaled', 'C_scaled', 'A_scaled:R_scaled', 'A_scaled:C_scaled', 'R_scaled:C_scaled', 'A_scaled:R_scaled:C_scaled']}

    for protocol, df_protocol in protocol_dict.items():

        fig = make_subplots(
            rows=1,
            cols=1,
            shared_xaxes=True,
        )

        df_protocol["sig"] = df_protocol["p.value"].apply(p_to_stars)
        df_protocol = df_protocol.query(f"term in {LMM_param_dict[protocol]}")

        y_max = df_protocol["estimate"].abs().max()

        fig.add_bar(
            x=df_protocol["term"],
            y=df_protocol["estimate"],
            text=df_protocol["sig"],
            textposition="outside",
            row=1,
            col=1
        )

        margin = y_max * 0.50   # 15% extra space

        fig.update_yaxes(
            range=[-y_max - margin, y_max + margin],
            row=1,
            col=1
        )

        fig.update_layout(
            title=f"ETmm {protocol}",
            template="simple_white",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            height=400,  # scale height
            width=500,
        )

        # fig.show()

        fig.write_html(
            os.path.join(path_export_res, f"{sujet}_{protocol}_LMM.html"),
            include_plotlyjs="cdn"
        )








################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        plot_median_ETmm_allcond(sujet)

        plot_REG_patient_wise(sujet)
        plot_STATS_REG_patient_wise(sujet)

    

        


    























