

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
######## EXPORT DF ########
################################

#ref_type = 'mn'
def export_Pxx_df(ref_type):
        
    #### load
    df_Pxx_allsujet = []

    for sujet in sujet_list:

        _df_sujet = []

        for cond in conditions:

            _xr = xr.open_dataarray(os.path.join(path_precompute, 'Pxx', f"{sujet}_{cond}_xr_Pxx_{ref_type}.nc"))
            _df_Pxx = _xr.to_dataframe().reset_index().dropna()
            _df_Pxx['cond'] = [cond] * _df_Pxx.shape[0]
            _df_sujet.append(_df_Pxx)

        _df_sujet = pd.concat(_df_sujet)

        respfeatures = get_respfeatures(sujet).rename(columns={'trial' : 'cycle_trial', 'cycle' : 'cycle_id'})

        _df_Pxx_ocratio = pd.merge(left=_df_sujet, right=respfeatures[['cycle_id', 'cycle_trial', 'oc_ratio']], how='inner', on=['cycle_id', 'cycle_trial'])

        df_Pxx_allsujet.append(_df_Pxx_ocratio)

    df_Pxx_allsujet = pd.concat(df_Pxx_allsujet)

    #### export
    df_Pxx_allsujet.to_excel(os.path.join(path_precompute, 'Pxx', 'df_R', f"df_R_Pxx.xlsx"))

    






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #ref_type = 'mn'
    for ref_type in ref_type_list:

        export_Pxx_df(ref_type)

