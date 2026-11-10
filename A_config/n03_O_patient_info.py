


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *


################################
######## TRIAL LIST ########
################################


trial_list_dict = {'P01' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P02' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P03' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P04' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P05' : ['oc300', 'oc600'],
                   'P06' : ['oc300'],
                   'P07' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P08' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P09' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P10' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P11' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P12' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P13' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P14' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P15' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   'P16' : ['oc300', 'oc600', 'oc300', 'oc600'],
                   }


shift_all_patient = {'P01' : {'trial0' : [-4.75], 'trial1' : [0], 'trial2' : [-4.9], 'trial3' : [0], },
                    'P02' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P03' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P04' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P05' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P06' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P07' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P08' : {'trial0' : [0], 'trial1' : [10.34], 'trial2' : [0], 'trial3' : [0], },
                    'P09' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P10' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P11' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                    'P12' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [-0.4], },
                    'P13' : {'trial0' : [-0.4], 'trial1' : [-0.5], 'trial2' : [-0.4], 'trial3' : [-0.5], },
                    'P14' : {'trial0' : [-0.675], 'trial1' : [-0.425], 'trial2' : [0], 'trial3' : [0], },
                    'P15' : {'trial0' : [0], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
                  }

# shift_all_patient = {'P01' : {'trial0' : [-4.75], 'trial1' : [0], 'trial2' : [-4.9], 'trial3' : [0], },
#                     'P02' : {'trial0' : [-17.51], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
#                     'P03' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
#                     'P04' : {'trial0' : [-8.55], 'trial1' : [0], 'trial2' : [0], 'trial3' : [0], },
#                           }


thresh_trial_start_stop_allpatient  = {'P01' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P02' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [], },
                    'P03' : {'trial0' : [0.5], 'trial1' : [0.25], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P04' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.75], },
                    'P05' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P06' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P07' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },     
                    'P08' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P09' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P10' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P11' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P12' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                    'P13' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.75], 'trial3' : [0.5], },
                    'P14' : {'trial0' : [0.5], 'trial1' : [0.5], 'trial2' : [0.5], 'trial3' : [0.5], },
                                              }



dict_chantrig_extraction_allpatient = {'P01' : ['CDIG_IN_1_Up', 'CDIG_IN_1_Down'], 'P02' : ['CDIG_IN_1_Up', 'CDIG_IN_1_Down'], 
                                       'P03' : 'CANALOG_IN_2', 'P04' : 'CANALOG_IN_1', 'P05' : 'CANALOG_IN_2',
                                       'P06' : 'CANALOG_IN_1', 'P07' : 'CANALOG_IN_1', 'P08' : 'CANALOG_IN_1', 'P09' : 'CANALOG_IN_1', 'P10' : 'CANALOG_IN_1',
                                       'P11' : 'CANALOG_IN_1', 'P12' : 'CANALOG_IN_1', 'P13' : 'CANALOG_IN_1', 'P14' : 'CANALOG_IN_1', 'P15' : 'CANALOG_IN_1',}




dict_chantrig_srate_extraction_allpatient = {'P01' : 'CDIG_IN_1_KHz', 'P02' : 'CDIG_IN_1_KHz', 
                                       'P03' : 'CANALOG_IN_2_KHz', 'P04' : 'CANALOG_IN_1_KHz', 'P05' : 'CANALOG_IN_2_KHz',
                                       'P06' : 'CANALOG_IN_1_KHz', 'P07' : 'CANALOG_IN_1_KHz', 'P08' : 'CANALOG_IN_1_KHz', 'P09' : 'CANALOG_IN_1_KHz', 'P10' : 'CANALOG_IN_1_KHz',
                                       'P11' : 'CANALOG_IN_1_KHz', 'P12' : 'CANALOG_IN_1_KHz', 'P13' : 'CANALOG_IN_1_KHz', 'P14' : 'CANALOG_IN_1_KHz', 'P15' : 'CANALOG_IN_1_KHz',}



correction_trig_SL_all_patient = {'P04' : np.array([[57900, 58100], [131300, 131350]]),
                               'P05' : np.array([[33700, 33800]]),
                               'P07' : np.array([[79000, 79500]]),
                               'P09' : np.array([[46000, 47000], [79800, 80000], [88000, 88200]]),
                               'P12' : np.array([[0, 9500]]),
                               'P13' : np.array([[2000, 4000], [44200, 44400], [13000, 16000], [32330, 32370], [38250, 38400], [102520, 102600], [168160, 168190]]),
                               'P14' : np.array([[70820, 70823]]),
                          }

correction_trig_iEEG_all_patient = {'P03' : {'trial0' : [], 'trial1' : [], 'trial2' : [], 'trial3' : [[467000,467004]]},
                                    'P07' : {'trial0' : [], 'trial1' : [], 'trial2' : [[153492,153500]], 'trial3' : []},
                                    'P09' : {'trial0' : [[474000,475500]], 'trial1' : [], 'trial2' : [], 'trial3' : []},
                                    'P10' : {'trial0' : [[458500,460540], [460550,461000]], 'trial1' : [], 'trial2' : [[167894,167900]], 'trial3' : [[71815, 71819], [596268, 596274], [656500, 657500], [730932, 730937]]},
                                    'P14' : {'trial0' : [], 'trial1' : [], 'trial2' : [[20000, 40000]], 'trial3' : []},
                          }


fasttrig_correction_allsujet = {'P01' : 2.5, 'P02' : 2.1, 'P03' : 2.5, 'P04' : 2.5, 'P05' : 1.5, 'P06' : 2.5, 'P07' : 1.75, 'P08' : 2,
                                'P09' : 2, 'P10' : 2, 'P11' : 2.5, 'P12' : 2, 'P13' : 2, 'P14' : 1.5, 'P15' : 2.5, 'P16' : 2.5,}



correction_cycle_detect_all_patient = { 'P02' : {'trial0' : [[19.9, 26.424], [81.46, 83.86]], 'trial1' : [[209.106, 210.137]], 'trial2' : [[247.46, 252.93]], 'trial3' : [[48.19, 55.43], [57.17, 61.24], [203.76, 204.9]]},
                                        'P03' : {'trial0' : [[277.88, 279.75]], 'trial1' : [], 'trial2' : [[205.66, 207.34]], 'trial3' : [[35.28, 38.19], [127.07, 129.18], [137.58, 139.69], [208.27, 209.60], [261.33, 263.14]]},
                                        'P07' : {'trial0' : [[14.26, 15.26], [34.03, 34.81], [66.376, 67.651], [105.913, 107.822], [142.149, 144.212], [150.89, 153.67], [173.16, 174.20], [289.89, 291.32]], 'trial1' : [[33.018, 35.546], [53.18, 54.75], [143.04, 146.35]], 'trial2' : [[56.49, 58.27], [171.872, 175.462], [241.21, 243.12]], 'trial3' : []},
                                        'P08' : {'trial0' : [[94.84, 105.92], [112.32, 113.626], [223.98, 225.799]], 'trial1' : [[166.38, 167.604], [212.23, 214.52]], 'trial2' : [[9.854, 11.891], [68.26, 71.27], [138.06, 142.05], [269.25, 273.887]], 'trial3' : [[128.93, 132.36], [145.631, 148.392]]},
                                        'P10' : {'trial0' : [], 'trial1' : [], 'trial2' : [[80.80, 86.81], [110.08, 113.80], [165.118, 169.300], [197.39, 201.96], [247.56, 251.09], [291.31, 296.44]], 'trial3' : [[36.56, 38.94], [47.83, 51.91], [141.36, 146.80], [264.55, 267.73]]},
                                        'P13' : {'trial0' : [[8.85, 10.511], [18.45, 20.67],  [256, 257]], 'trial1' : [], 'trial2' : [[44, 45.116], [72, 74]], 'trial3' : []},
                                        'P14' : {'trial0' : [[92.89, 94.87]], 'trial1' : [], 'trial2' : [], 'trial3' : []},
                          }



dict_nchan_allsujet = {'P01' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_01___Central', 'CMacro_LFP_02___Posteriomedial', 'CMacro_LFP_03___Anteriomedial', 'CMacro_LFP_04___Posteriolateral', 'CMacro_LFP_05___Anteriolateral']}, 
                       'P02' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_AP_T', 'CMacro_LFP_NP1', 'CMacro_LFP_NP2']}, 
                       'P03' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_01___Central', 'CMacro_LFP_02___Posteriomedial', 'CMacro_LFP_03___Anteriomedial', 'CMacro_LFP_04___Anteriolateral', 'CMacro_LFP_05___Posteriolateral']}, 
                       'P04' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_01___Central']}, 
                       'P05' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_AP_T', 'CMacro_LFP_NP1', 'CMacro_LFP_NP2']}, 
                       'P07' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_01___Central', 'CMacro_LFP_02___Posteriomedial', 'CMacro_LFP_03___Anteriomedial', 'CMacro_LFP_04___Anteriolateral', 'CMacro_LFP_05___Posteriolateral']}, 
                       'P08' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_01___Central']},
                       'P09' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_01___Central', 'CMacro_LFP_02___Posteriomedial', 'CMacro_LFP_03___Anteriomedial', 'CMacro_LFP_04___Anteriolateral', 'CMacro_LFP_05___Posteriolateral']}, 
                       'P10' : {'ECOG' : np.concat([np.arange(start=1, stop=7, step=1), np.arange(start=9, stop=15, step=1)]), 'MACRO' : ['CMacro_LFP_01___Central', 'CMacro_LFP_02___Posteriomedial', 'CMacro_LFP_03___Anteriomedial', 'CMacro_LFP_04___Anteriolateral', 'CMacro_LFP_05___Posteriolateral']}, 
                       'P12' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_AP_T', 'CMacro_LFP_NP1', 'CMacro_LFP_NP2']}, 
                       'P13' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_AP_T', 'CMacro_LFP_NP1', 'CMacro_LFP_NP2']},
                       'P14' : {'ECOG' : np.arange(start=1, stop=17, step=1), 'MACRO' : ['CMacro_LFP_AP_T', 'CMacro_LFP_NP1', 'CMacro_LFP_NP2']},
                       'P15' : {'ECOG' : 12, 'MACRO' : 1},
                       }



# no zero based
chan_ref_allsujet = {'P01' : 1, 'P02' : 1, 'P03' : 1, 'P04' : 1, 'P05' : 9, 'P06' : None, 'P07' : 1, 'P08' : 'subG',}



