from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import os
import numpy as np
import xarray as xr
import scipy
import matplotlib.pyplot as plt
from fractions import Fraction


debug = False


################################
######## PLOT FUNCTION ########
################################

def plot_xr_channels_matplotlib(xr_data, chan_start=0, n_chan=8, target_srate=100, srate=500, offset=4, include_resp=True,
                                norm_method_ieeg=None, norm_method_resp="zscore", resp_amp=20, vlines=None, title=None, figsize=(16, 10),):

    #### channel names
    chan_names = xr_data["chan"].values.astype(str)

    #### respiration
    resp_exists = "resp" in chan_names
    data_chan_names = chan_names[chan_names != "resp"]

    #### select channels
    chan_stop = min(chan_start + n_chan, len(data_chan_names))
    chan_sel = data_chan_names[chan_start:chan_stop]

    ################################
    ######## RESAMPLING ########
    ################################

    ratio = Fraction(target_srate / srate).limit_denominator(1000)
    up = ratio.numerator
    down = ratio.denominator

    ################################
    ######## NORMALIZATION ########
    ################################

    def normalize_signal(x, norm_method):

        if norm_method == "zscore":
            return scipy.stats.zscore(x)

        elif norm_method == "robust":
            med = np.median(x)
            mad = scipy.stats.median_abs_deviation(x)

            if mad == 0:
                return np.full_like(x, np.nan, dtype=float)

            return (x - med) * 0.6745 / mad

        elif norm_method is None:
            return x

        else:
            raise ValueError("norm_method must be 'zscore', 'robust', or None")

    ################################
    ######## FIGURE ########
    ################################

    fig, ax = plt.subplots(figsize=figsize)

    trace_i = 0
    yticks = []
    yticklabels = []
    time_down = None

    ################################
    ######## RESPIRATION ########
    ################################

    if include_resp and resp_exists:

        resp = xr_data.sel(chan="resp").values

        if down != 1 or up != 1:
            resp_down = scipy.signal.resample_poly(resp, up=up, down=down)
        else:
            resp_down = resp

        if time_down is None:
            time_down = np.arange(resp_down.size) / target_srate

        y = normalize_signal(resp_down, norm_method_resp) * resp_amp
        y_offset = y - trace_i * offset

        ax.plot(time_down, y_offset, linewidth=1.2)

        yticks.append(-trace_i * offset)
        yticklabels.append("resp")

        trace_i += 1

    ################################
    ######## IEEG CHANNELS ########
    ################################

    for chan in chan_sel:

        x = xr_data.sel(chan=chan).values

        if down != 1 or up != 1:
            x_down = scipy.signal.resample_poly(x, up=up, down=down)
        else:
            x_down = x

        if time_down is None:
            time_down = np.arange(x_down.size) / target_srate

        y = normalize_signal(x_down, norm_method_ieeg)
        y_offset = y - trace_i * offset

        ax.plot(time_down, y_offset, linewidth=0.8)

        yticks.append(-trace_i * offset)
        yticklabels.append(chan)

        trace_i += 1

    ################################
    ######## VERTICAL LINES ########
    ################################

    if vlines is not None:
        for line in vlines:
            ax.axvline(x=line, color="red", linestyle="--", linewidth=1.2)

    ################################
    ######## LAYOUT ########
    ################################

    if title is None:
        title = f"Channels {chan_start}–{chan_stop - 1}"

    ax.set_title(title, fontsize=16, fontweight="bold")
    ax.set_xlabel("Time (s)", fontsize=14)
    ax.set_ylabel("Channels", fontsize=14)

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=10)

    ax.tick_params(axis="x", labelsize=11)

    if time_down is not None:
        ax.set_xlim(time_down[0], time_down[-1])

    ax.grid(axis="x", alpha=0.2)

    plt.tight_layout()
    plt.show()

    return fig, ax













################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    ################################
    ######## CONFIG ########
    ################################

    
    ref_type = 'mn'

    #sujet_list
    sujet = 'P02'
    trial_i = 0

    chanlist = get_chanlist_raw(sujet, ref_type)

    srate = 500
    target_srate = 100
    
    offset = 800

    ################################
    ######## LOAD DATA ########
    ################################

    xr_data = xr.open_dataarray(os.path.join(path_prep, 'export', sujet, f"{sujet}_trial{trial_i+1}_raw_{ref_type}.nc"))


    # xr_data_preproc = xr.open_dataarray(os.path.join(path_prep, 'export', sujet, f"{sujet}_trial{trial_i+1}_preproc.nc"))
    # xr_data_bi = xr.open_dataarray(os.path.join(path_prep, 'export', sujet, f"{sujet}_trial{trial_i+1}_raw_bi.nc"))
    # xr_data_preproc_bi = xr.open_dataarray(os.path.join(path_prep, 'export', sujet, f"{sujet}_trial{trial_i+1}_preproc_bi.nc"))



    ################################
    ######## VISU ########
    ################################

    plot_xr_channels_matplotlib(xr_data, chan_start=0, n_chan=len(chanlist), target_srate=target_srate, srate=srate, offset=offset, include_resp=True,
                                norm_method_ieeg=None, norm_method_resp=None, resp_amp=40, vlines=None, title=f"{sujet} trial {trial_i+1} RAW")

    

    # plot_xr_channels_matplotlib(xr_data_preproc, chan_start=0, n_chan=16, target_srate=target_srate, srate=srate, offset=offset, include_resp=True,
    #                             norm_method_ieeg=None, norm_method_resp="zscore", resp_amp=40, vlines=None, title=f"{sujet} trial {trial_i+1} PREPROC")

    # plot_xr_channels_matplotlib(xr_data_bi, chan_start=0, n_chan=16, target_srate=target_srate, srate=srate, offset=offset, include_resp=True,
    #                                 norm_method_ieeg=None, norm_method_resp="zscore", resp_amp=40, vlines=None, title=f"{sujet} trial {trial_i+1} BI")
    
    # plot_xr_channels_matplotlib(xr_data_preproc_bi, chan_start=0, n_chan=16, target_srate=target_srate, srate=srate, offset=offset, include_resp=True,
    #                             norm_method_ieeg=None, norm_method_resp="zscore", resp_amp=40, vlines=None, title=f"{sujet} trial {trial_i+1} PREPROC BI")



    ########################################
    ######## COMPARE ONE CHANNEL ########
    ########################################

    chan_i = 2

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.plot(xr_data.data[chan_i], label="raw", alpha=0.8)
    ax.plot(xr_data_preproc.data[chan_i], label="preproc", alpha=0.8)
    ax.plot(xr_data_bi.data[chan_i], label="raw_bi", alpha=0.8)
    ax.plot(xr_data_preproc_bi.data[chan_i], label="preproc_bi", alpha=0.8)    

    ax.legend()
    ax.set_title(f"{sujet} channel comparison")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Amplitude")

    plt.tight_layout()
    plt.show()


    # ################################
    # ######## AVERAGE REFERENCE ########
    # ################################

    # chanlist, localist = get_chanlocalist(sujet, 'mn')

    # xr_ieeg = xr_data.sel(chan=xr_data["chan"] != "resp")
    # av_ref = xr_ieeg.median(dim="chan").data

    # av_ref_data = []

    # for chan in chanlist:
    #     x = xr_data.sel(chan=chan).data
    #     x_ref = x - av_ref
    #     av_ref_data.append(x_ref)

    # av_ref_data = np.vstack(av_ref_data)

    # xr_av_ref = xr.DataArray(
    #     av_ref_data,
    #     dims=("chan", "time"),
    #     coords={"chan": chanlist, "time": xr_data["time"].values}
    # )

    # plot_xr_channels_matplotlib(
    #     xr_av_ref,
    #     chan_start=0,
    #     n_chan=16,
    #     target_srate=target_srate,
    #     srate=srate,
    #     offset=offset,
    #     include_resp=False,
    #     norm_method_ieeg=None,
    #     norm_method_resp="zscore",
    #     vlines=None,
    #     title=f"{sujet} trial {trial_i+1} AVERAGE REF"
    # )


    