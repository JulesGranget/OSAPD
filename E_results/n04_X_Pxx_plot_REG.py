



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




def plot_REG_patient_wise(sujet):

    #### config
    path_export_res = os.path.join(path_results, 'Pxx', 'REG', 'patientwise')

    #### load df
    path_load_reg = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx')

    df_reg_VCCP = [pd.read_excel(os.path.join(path_load_reg, f"{sujet}_df_reg_Pxx_{band}_VCCP_R.xlsx")) for band in freq_band_dict]
    df_reg_VCCP = pd.concat(df_reg_VCCP).drop(columns=['Unnamed: 0'])
    df_reg_VCCP = df_reg_VCCP.query(f"ROI in {ROI_short_list}")

    #### plot
    for ROI in ROI_short_list:

        df_plot = df_reg_VCCP.query(f"ROI == '{ROI}'")
        col_melt = ['A_scaled', 'R_scaled', 'C_scaled']
        df_plot = pd.melt(df_plot, id_vars=[col for col in df_plot.columns if col not in col_melt], value_vars=col_melt, var_name='factor', value_name='rscore')

        fig = sns.lmplot(data=df_plot, x='rscore', y='Pxx', hue='phase_cycle', row='band', col='factor', sharey=False, sharex=False)
        fig.fig.subplots_adjust(top=0.92, wspace=0.35, hspace=0.45)
        fig.fig.suptitle(f"VCCP {ROI}")
        
        # fig.fig.show()

        fig.fig.savefig(os.path.join(path_export_res, f"{sujet}_VCCP_REG_{ROI}.png"))









################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        plot_REG_patient_wise(sujet)

    

        


    























