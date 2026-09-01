

import joblib


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import joblib
import statsmodels.api as sm
from io import StringIO
from statsmodels.tools.sm_exceptions import ValueWarning
import warnings





################################
######## EXPORT DF COND ########
################################



#sujet = sujet_list[0]
def export_Pxx_cond_df_sujet(sujet):

    for band in freq_band_dict:

        print(band)
        
        #### load
        path_load_Pxx = os.path.join(path_precompute, 'Pxx')

        _xr = xr.open_dataarray(os.path.join(path_load_Pxx, f"{sujet}_xr_Pxx.nc")).sel(band=band).rename({'phase' : 'phase_cycle'})

        #### generate df
        df_band_Pxx_sujet = _xr.to_dataframe().reset_index()
        
        #### clean
        nan_exclude_vec = ~df_band_Pxx_sujet['Pxx'].isna()
        df_band_Pxx_sujet = df_band_Pxx_sujet[nan_exclude_vec]
        df_band_Pxx_sujet = df_band_Pxx_sujet.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

        #### prepare cond
        df_band_Pxx_sujet_attention = df_band_Pxx_sujet.query(f"cond in ['RB', 'AoRoCo']")
        df_band_Pxx_sujet_HV = df_band_Pxx_sujet.query(f"cond in ['RB', 'HV']")
        df_band_Pxx_sujet_VCCP = df_band_Pxx_sujet.query(f"cond not in ['RB', 'HV']")

        A_cond, R_cond, C_cond = [], [], []
        for row_i, row_val in df_band_Pxx_sujet_VCCP.iterrows():

            A_cond.append(row_val['cond'][1])
            R_cond.append(row_val['cond'][3])
            C_cond.append(row_val['cond'][5])

        df_band_Pxx_sujet_VCCP['A'], df_band_Pxx_sujet_VCCP['R'], df_band_Pxx_sujet_VCCP['C'] = A_cond, R_cond, C_cond
        df_band_Pxx_sujet_VCCP = df_band_Pxx_sujet_VCCP.drop(columns=['cond'])

        df_band_Pxx_sujet_attention['sujet'] = [sujet] * df_band_Pxx_sujet_attention.shape[0]
        df_band_Pxx_sujet_HV['sujet'] = [sujet] * df_band_Pxx_sujet_HV.shape[0]
        df_band_Pxx_sujet_VCCP['sujet'] = [sujet] * df_band_Pxx_sujet_VCCP.shape[0]

        #### calibrate for different protocols
        df_band_Pxx_sujet_VCCP['C'] = df_band_Pxx_sujet_VCCP['C'].replace({'c' : 'o'})
        df_band_Pxx_sujet_attention['cond'] = df_band_Pxx_sujet_attention['cond'].replace({'AoRoCo' : 'CB'})

        #### export
        path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx')

        df_band_Pxx_sujet_attention.to_excel(os.path.join(path_export_df, f"{sujet}_df_Pxx_{band}_attention_R.xlsx"))
        df_band_Pxx_sujet_HV.to_excel(os.path.join(path_export_df, f"{sujet}_df_Pxx_{band}_HV_R.xlsx"))
        df_band_Pxx_sujet_VCCP.to_excel(os.path.join(path_export_df, f"{sujet}_df_Pxx_{band}_VCCP_R.xlsx"))






#sujet = sujet_list[0]
def export_Cxy_cond_df_sujet(sujet):
        
    #### load
    path_load_Cxy = os.path.join(path_precompute, 'Cxy')

    chanlist, _, _, localist = get_chanlist_localist(sujet)
    _xr_Cxy = xr.open_dataarray(os.path.join(path_load_Cxy, f"{sujet}_Cxy.nc"))

    ### extract Cxy around respi
    respfeature = get_respfeatures_main(sujet)

    trial_list = get_alltrials_sujet(sujet)

    resp_med_allcond = {}

    #cond = conditions[-1]
    for cond in conditions:

        trial_list_cond = [trial for trial in trial_list if trial.find(cond) != -1]

        _resp_med_cond = []

        for trial_i, trial in enumerate(trial_list_cond):

            _respfeature_trial = respfeature.query(f"cond == '{cond}' and trial == {trial_i+1}")
            _resp_med_cond.append(_respfeature_trial['cycle_freq'].median())

        resp_med_allcond[cond] = np.median(_resp_med_cond)

    data_Cxy_med = []
    hzCxy = _xr_Cxy['freq'].values

    for cond in conditions:

        sel_vec_hzCxy = (hzCxy < (resp_med_allcond[cond] + Cxy_extraction_range)) & ((resp_med_allcond[cond] - Cxy_extraction_range) < hzCxy) 
        data_Cxy_med.append(_xr_Cxy.sel(cond=cond, freq=sel_vec_hzCxy).median('freq').values)

    data_Cxy_med = np.stack(data_Cxy_med)

    xr_coords = {'cond' : conditions, 'chan' : chanlist, 'ROI' : ('chan', localist['loca'].values)}
    _xr_Cxy_med = xr.DataArray(data_Cxy_med, dims=['cond', 'chan'], coords=xr_coords) 

    #### generate df
    df_Cxy_sujet_raw = _xr_Cxy_med.to_dataframe(name='Cxy').reset_index()
    
    #### clean
    df_Cxy_sujet = df_Cxy_sujet_raw.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

    #### prepare cond
    df_Cxy_sujet_attention = df_Cxy_sujet.query(f"cond in ['RB', 'AoRoCo']")
    df_Cxy_sujet_HV = df_Cxy_sujet.query(f"cond in ['RB', 'HV']")
    df_Cxy_sujet_protocol = df_Cxy_sujet.query(f"cond not in ['RB', 'HV']")

    A_cond, R_cond, C_cond = [], [], []
    for row_i, row_val in df_Cxy_sujet_protocol.iterrows():

        A_cond.append(row_val['cond'][1])
        R_cond.append(row_val['cond'][3])
        C_cond.append(row_val['cond'][5])

    df_Cxy_sujet_protocol['A'], df_Cxy_sujet_protocol['R'], df_Cxy_sujet_protocol['C'] = A_cond, R_cond, C_cond
    df_Cxy_sujet_protocol = df_Cxy_sujet_protocol.drop(columns=['cond'])

    df_Cxy_sujet_attention['sujet'] = [sujet] * df_Cxy_sujet_attention.shape[0]
    df_Cxy_sujet_HV['sujet'] = [sujet] * df_Cxy_sujet_HV.shape[0]
    df_Cxy_sujet_protocol['sujet'] = [sujet] * df_Cxy_sujet_protocol.shape[0]

    #### calibrate for different protocols
    df_Cxy_sujet_protocol['C'] = df_Cxy_sujet_protocol['C'].replace({'c' : 'o'})
    df_Cxy_sujet_attention['cond'] = df_Cxy_sujet_attention['cond'].replace({'AoRoCo' : 'CB'})

    #### export
    path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Cxy')

    df_Cxy_sujet_attention.to_excel(os.path.join(path_export_df, f"{sujet}_df_Cxy_attention_R.xlsx"))
    df_Cxy_sujet_HV.to_excel(os.path.join(path_export_df, f"{sujet}_df_Cxy_HV_R.xlsx"))
    df_Cxy_sujet_protocol.to_excel(os.path.join(path_export_df, f"{sujet}_df_Cxy_protocol_R.xlsx"))









################################
######## EXPORT DF REG ########
################################






#sujet = sujet_list[0]
def export_Pxx_reg_df_sujet(sujet):

    #### config
    sel_vec_rf_metric_scaled = ['sujet', 'cond', 'cycle', 'A', 'R', 'C', 'A_scaled', 'R_scaled', 'C_scaled']

    #band = 'theta'
    for band in freq_band_dict:

        print(band)
        
        #### load
        path_load_Pxx = os.path.join(path_precompute, 'Pxx')
        path_load_rf_scaled = os.path.join(path_precompute, 'RESPI', 'respfeatures', sujet)
        
        _xr = xr.open_dataarray(os.path.join(path_load_Pxx, f"{sujet}_xr_Pxx.nc")).sel(band=band).rename({'phase' : 'phase_cycle'})

        #### generate df
        df_band_Pxx_sujet = _xr.to_dataframe().reset_index()
        
        #### clean
        nan_exclude_vec = ~df_band_Pxx_sujet['Pxx'].isna()
        df_band_Pxx_sujet = df_band_Pxx_sujet[nan_exclude_vec]
        df_band_Pxx_sujet = df_band_Pxx_sujet.query(f"ROI in {ROI_short_list}").reset_index(drop=True)

        #### prepare cond
        df_band_Pxx_sujet_attention = df_band_Pxx_sujet.query(f"cond in ['RB', 'AoRoCo']")
        df_band_Pxx_sujet_attention['cond'] = df_band_Pxx_sujet_attention['cond'].replace({'AoRoCo' : 'CB'})
        df_band_Pxx_sujet_VCCP = df_band_Pxx_sujet.query(f"cond not in ['RB', 'HV']")

        #### load Pxx and rf scaled
        rf_SD_scaled_VCCP = pd.read_excel(os.path.join(path_load_rf_scaled, f"respfeatures_scaled_VCCP.xlsx")).query(f"cond not in ['RB', 'HV', 'excluded']")[sel_vec_rf_metric_scaled]
        rf_SD_scaled_RB = pd.read_excel(os.path.join(path_load_rf_scaled, f"respfeatures_scaled_RB.xlsx"))

        #### prepare attention
        rf_SD_scaled_attention = rf_SD_scaled_RB.query(f"cond in ['RB', 'AoRoCo']")
        rf_SD_scaled_attention['cond'] = rf_SD_scaled_attention['cond'].replace({'AoRoCo' : 'CB'})
        rf_SD_scaled_attention = rf_SD_scaled_attention[sel_vec_rf_metric_scaled]

        #### merge
        df_reg_VCCP = pd.merge(df_band_Pxx_sujet_VCCP, rf_SD_scaled_VCCP, on=['sujet', 'cond', 'cycle'], how='left')
        df_reg_attention = pd.merge(df_band_Pxx_sujet_attention, rf_SD_scaled_attention, on=['sujet', 'cond', 'cycle'], how='left')

        #### calibrate for different protocols
        df_reg_VCCP['C'] = df_reg_VCCP['C'].replace({'c' : 'o'})        

        #### export
        path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'Pxx')

        df_reg_attention.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_Pxx_{band}_attention_R.xlsx"))
        df_reg_VCCP.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_Pxx_{band}_VCCP_R.xlsx"))








#sujet = sujet_list[0]
def export_ET_reg_df_sujet(sujet):

    #### config
    sel_vec_rf_metric_scaled = ['sujet', 'cond', 'cycle', 'A_scaled', 'R_scaled', 'C_scaled']
    sel_vec_ETmm_metric_scaled = ['sujet', 'cond', 'cycle', 'ETmm_scaled']

    #### load
    path_load_ETmm = os.path.join(path_precompute, 'ET')
    path_load_rf = os.path.join(path_precompute, 'RESPI', 'respfeatures', sujet)
    
    #### merge
    df_scaled_VCCP_ETmm = pd.read_excel(os.path.join(path_load_ETmm, f"df_rf_ETmm_scaled_VCCP.xlsx")).drop(columns=['Unnamed: 0']).query(f"cond not in ['RB', 'HV']")[sel_vec_ETmm_metric_scaled]
    df_scaled_VCCP_rf = pd.read_excel(os.path.join(path_load_rf, f"respfeatures_scaled_VCCP.xlsx")).drop(columns=['Unnamed: 0']).query(f"cond not in ['RB', 'HV']")[sel_vec_rf_metric_scaled]
    df_scaled_VCCP = pd.merge(df_scaled_VCCP_ETmm, df_scaled_VCCP_rf, how='inner', on=['sujet', 'cond', 'cycle'])

    df_scaled_HV_ETmm = pd.read_excel(os.path.join(path_load_ETmm, f"df_rf_ETmm_scaled_RB.xlsx")).drop(columns=['Unnamed: 0']).query(f"cond in ['RB', 'HV']")[sel_vec_ETmm_metric_scaled]
    df_scaled_HV_rf = pd.read_excel(os.path.join(path_load_rf, f"respfeatures_scaled_RB.xlsx")).drop(columns=['Unnamed: 0']).query(f"cond in ['RB', 'HV']")[sel_vec_rf_metric_scaled]
    df_scaled_HV = pd.merge(df_scaled_HV_ETmm, df_scaled_HV_rf, how='inner', on=['sujet', 'cond', 'cycle'])

    df_scaled_CB_ETmm = pd.read_excel(os.path.join(path_load_ETmm, f"df_rf_ETmm_scaled_RB.xlsx")).drop(columns=['Unnamed: 0']).query(f"cond in ['RB', 'AoRoCo']")[sel_vec_ETmm_metric_scaled]
    df_scaled_CB_rf = pd.read_excel(os.path.join(path_load_rf, f"respfeatures_scaled_RB.xlsx")).drop(columns=['Unnamed: 0']).query(f"cond in ['RB', 'AoRoCo']")[sel_vec_rf_metric_scaled]
    df_scaled_CB = pd.merge(df_scaled_CB_ETmm, df_scaled_CB_rf, how='inner', on=['sujet', 'cond', 'cycle'])
    df_scaled_CB['cond'] = df_scaled_CB['cond'].replace({'AoRoCo' : 'CB'})    

    #### export
    path_export_df = os.path.join(path_precompute, 'EXPORT_DF', 'patient_wise', 'ETmm')

    df_scaled_VCCP.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_ETmm_VCCP_R.xlsx"))
    df_scaled_HV.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_ETmm_HV_R.xlsx"))
    df_scaled_CB.to_excel(os.path.join(path_export_df, f"{sujet}_df_reg_ETmm_CB_R.xlsx"))

















################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        export_Pxx_cond_df_sujet(sujet)
        export_Cxy_cond_df_sujet(sujet)
        
        export_Pxx_reg_df_sujet(sujet)
        export_ET_reg_df_sujet(sujet)







                        