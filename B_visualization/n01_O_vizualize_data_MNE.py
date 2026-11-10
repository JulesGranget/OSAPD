from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import os
import numpy as np
import xarray as xr
import scipy
import matplotlib.pyplot as plt
import mne


debug = False


################################
######## MNE PLOT FUNCTION #####
################################

def plot_xr_channels_mne(xr_data, srate=500, target_srate=100, chan_start=0, n_chan=None, include_resp=True,
                         seeg_scaling=800, resp_scaling=40, duration=20, title=None):

    #### channel names
    chan_names = xr_data["chan"].values.astype(str).tolist()

    ################################
    ######## REMOVE RESP ###########
    ################################

    if not include_resp and "resp" in chan_names:
        xr_data = xr_data.sel(chan=xr_data["chan"] != "resp")
        chan_names = xr_data["chan"].values.astype(str).tolist()

    ################################
    ######## CHANNEL TYPES #########
    ################################

    ch_types = ["resp" if chan == "resp" else "seeg" for chan in chan_names]

    ################################
    ######## CREATE MNE RAW ########
    ################################

    info = mne.create_info(ch_names=chan_names, sfreq=srate, ch_types=ch_types)

    raw = mne.io.RawArray(
        xr_data.values,
        info,
        verbose=False
    )

    ################################
    ######## RESAMPLING ############
    ################################

    if target_srate != srate:
        raw.resample(target_srate)

    ################################
    ######## SELECT CHANNELS #######
    ################################

    data_chan_names = [chan for chan in raw.ch_names if chan != "resp"]

    if n_chan is None:
        n_chan = len(data_chan_names)

    chan_stop = min(chan_start + n_chan, len(data_chan_names))
    chan_sel = data_chan_names[chan_start:chan_stop]

    if include_resp and "resp" in raw.ch_names:
        picks = ["resp"] + chan_sel
    else:
        picks = chan_sel

    ################################
    ######## TITLE #################
    ################################

    if title is None:
        title = f"Channels {chan_start}–{chan_stop - 1}"

    ################################
    ######## MNE BROWSER ###########
    ################################

    raw.plot(
        picks=picks,
        duration=duration,
        n_channels=len(picks),
        scalings={"seeg": seeg_scaling, "resp": resp_scaling},
        remove_dc=False,
        group_by="original",
        title=title,
        show_scrollbars=True,
        show_scalebars=True,
        block=True
    )

    return raw













################################
######## EXECUTE ###############
################################

if __name__ == '__main__':

    ################################
    ######## CONFIG ###############
    ################################

    ref_type = 'mn'

    #sujet_list
    sujet = 'P14'
    trial_i = 0

    chanlist = get_chanlist_raw(sujet, ref_type)

    srate = 500
    target_srate = 100

    seeg_scaling = 800
    resp_scaling = 40

    ################################
    ######## LOAD DATA ############
    ################################

    xr_data = xr.open_dataarray(os.path.join(path_prep, 'export', sujet, f"{sujet}_trial{trial_i+1}_raw_{ref_type}.nc"))


    ################################
    ######## VISU RAW #############
    ################################

    raw_mne = plot_xr_channels_mne(xr_data, srate=srate, target_srate=target_srate, chan_start=0, n_chan=len(chanlist), include_resp=True,
        seeg_scaling=seeg_scaling, resp_scaling=resp_scaling, duration=20, title=f"{sujet} trial{trial_i} RAW")










    ################################
    ######## VISU PREPROC #########
    ################################

    # raw_preproc_mne = plot_xr_channels_mne(
    #     xr_data_preproc,
    #     srate=srate,
    #     target_srate=target_srate,
    #     chan_start=0,
    #     n_chan=len(chanlist),
    #     include_resp=True,
    #     seeg_scaling=seeg_scaling,
    #     resp_scaling=resp_scaling,
    #     duration=20,
    #     title=f"{sujet} trial {trial_i+1} PREPROC"
    # )



    ################################
    ######## VISU BIPOLAR #########
    ################################

    # raw_bi_mne = plot_xr_channels_mne(
    #     xr_data_bi,
    #     srate=srate,
    #     target_srate=target_srate,
    #     chan_start=0,
    #     n_chan=xr_data_bi.sizes["chan"],
    #     include_resp=True,
    #     seeg_scaling=seeg_scaling,
    #     resp_scaling=resp_scaling,
    #     duration=20,
    #     title=f"{sujet} trial {trial_i+1} BI"
    # )



    ################################
    ######## VISU PREPROC BI ######
    ################################

    # raw_preproc_bi_mne = plot_xr_channels_mne(
    #     xr_data_preproc_bi,
    #     srate=srate,
    #     target_srate=target_srate,
    #     chan_start=0,
    #     n_chan=xr_data_preproc_bi.sizes["chan"],
    #     include_resp=True,
    #     seeg_scaling=seeg_scaling,
    #     resp_scaling=resp_scaling,
    #     duration=20,
    #     title=f"{sujet} trial {trial_i+1} PREPROC BI"
    # )



    ########################################
    ######## COMPARE ONE CHANNEL ##########
    ########################################

    # chan_i = 2

    # fig, ax = plt.subplots(figsize=(14, 6))

    # ax.plot(xr_data.data[chan_i], label="raw", alpha=0.8)
    # ax.plot(xr_data_preproc.data[chan_i], label="preproc", alpha=0.8)
    # ax.plot(xr_data_bi.data[chan_i], label="raw_bi", alpha=0.8)
    # ax.plot(xr_data_preproc_bi.data[chan_i], label="preproc_bi", alpha=0.8)

    # ax.legend()
    # ax.set_title(f"{sujet} channel comparison")
    # ax.set_xlabel("Samples")
    # ax.set_ylabel("Amplitude")

    # plt.tight_layout()
    # plt.show()