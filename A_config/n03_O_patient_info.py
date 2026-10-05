


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *


################################
######## TRIAL LIST ########
################################


trial_list_dict = {'P01' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P02' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P03' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P04' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   }


shift_all_patient = {'P01' : {'trial0' : [-4.75], 'trial1' : [0], 'trial2' : [-4.9], 'trial3' : [0], },
                    'P02' : {'trial0' : [-17.51], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
                    'P03' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
                    'P04' : {'trial0' : [-8.55], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                          }


thresh_trial_start_stop_allpatient  = {'P01' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
                    'P02' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
                    'P03' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
                    'P04' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.75], },
                          }

correction_trig_all_patient = {'P04' : np.array([[57900, 58100], [131300, 131350]]),
                          }

# no zero based
chan_ref_allsujet = {'P01' : 1, 'P02' : 1, 'P03' : 1, 'P04' : 1,}



