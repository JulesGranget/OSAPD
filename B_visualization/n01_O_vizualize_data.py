
from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import plotly.graph_objects as go


debug = False








################################
######## INSPECT DATA ########
################################



def plot_xr_channels(xr_data, chan_start=0, n_chan=8, target_srate=100, srate=500, offset=4, include_resp=True, norm_method_ieeg=None, norm_method_resp="zscore", vlines=None, title=None,):


    #### channel names
    chan_names = xr_data["chan"].values.astype(str)

    #### separate respiration
    resp_exists = "resp" in chan_names

    data_chan_names = chan_names[chan_names != "resp"]

    #### select channel subset
    chan_stop = min(
        chan_start + n_chan,
        len(data_chan_names)
    )

    chan_sel = data_chan_names[
        chan_start:chan_stop
    ]

    #### downsampling ratio
    from fractions import Fraction

    ratio = Fraction(
        target_srate / srate
    ).limit_denominator(1000)

    up = ratio.numerator
    down = ratio.denominator

    time_down = None


    #### normalization function
    def normalize_signal(x, norm_method):

        if norm_method == "zscore":
            return scipy.stats.zscore(x)

        elif norm_method == "robust":

            med = np.median(x)

            mad = scipy.stats.median_abs_deviation(x)

            if mad == 0:
                return x * np.nan

            return (
                (x - med)
                * 0.6745
                / mad
            )

        elif norm_method is None:
            return x

        else:
            raise ValueError(
                "norm must be 'zscore', 'robust', or None"
            )

    #### figure
    fig = go.Figure()

    trace_i = 0

    ################################
    ######## RESPIRATION ########
    ################################

    if include_resp and resp_exists:

        resp = xr_data.sel(chan="resp").values

        if down != 1 or up != 1:
            resp_down = scipy.signal.resample_poly(
                resp,
                up=up,
                down=down
            )
        else:
            resp_down = resp

        if time_down is None:
            time_down = (
                np.arange(resp_down.size)
                / target_srate
            )

        y = normalize_signal(resp_down, norm_method_resp)

        fig.add_trace(
            go.Scattergl(
                x=time_down,
                y=y,
                mode="lines",
                name="resp"
            )
        )

        trace_i += 1

    ################################
    ######## IEEG CHANNELS ########
    ################################

    for chan_i, chan in enumerate(chan_sel):

        x = xr_data.sel(chan=chan).values

        if down != 1 or up != 1:
            x_down = scipy.signal.resample_poly(
                x,
                up=up,
                down=down
            )
        else:
            x_down = x

        if time_down is None:
            time_down = (
                np.arange(x_down.size)
                / target_srate
            )

        y = normalize_signal(x_down, norm_method_ieeg)

        y = y - trace_i * offset

        fig.add_trace(
            go.Scattergl(
                x=time_down,
                y=y,
                mode="lines",
                name=chan
            )
        )

        trace_i += 1

    ################################
    ######## VERTICAL LINES ########
    ################################

    if vlines is not None:

        for line in vlines:

            fig.add_vline(
                x=line,
                line_width=1.5,
                line_dash="dash",
                line_color="red"
            )

    ################################
    ######## LAYOUT ########
    ################################

    if title is None:
        title = (
            f"Channels {chan_start}"
            f"–{chan_stop - 1}"
        )

    fig.update_layout(
        template="simple_white",
        title=title,
        xaxis_title="Time (s)",
        yaxis_title="Amplitude (stacked)",
        height=max(
            500,
            100 * (len(chan_sel) + int(include_resp and resp_exists))
        ),
        hovermode="x unified"
    )

    fig.show()

    return fig












################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    #### load
    sujet = 'P01'
    trial_i = 0
    ref_type = 'mn'

    chanlist, localist = get_chanlocalist(sujet, ref_type)

    xr_data = xr.open_dataarray(os.path.join(path_prep, 'export', sujet, f"{sujet}_trial{trial_i+1}_raw.nc"))

    #### raw
    plot_xr_channels(xr_data, chan_start=0, n_chan=16, target_srate=100, srate=500, offset=4, include_resp=True, norm_method_ieeg=None, norm_method_resp="zscore", vlines=None, title=None,)

    #### average reref
    av_ref = xr_data.sel(chan=xr_data['chan'] != 'resp').median(dim='chan').data

    av_ref_data = [xr_data.sel(chan=chan).data - av_ref for chan in chanlist]
    av_ref_data = np.vstack(av_ref_data)

    xr_av_ref = xr.DataArray(
        av_ref_data,
        dims=("chan", "time"),
        coords={
            "chan": chanlist,
            "time": xr_data["time"].values
        }
    )

    plot_xr_channels(xr_av_ref, chan_start=0, n_chan=6, target_srate=100, srate=500, offset=4, include_resp=True, norm_resp="zscore", vlines=None, title=None,)

    #### biploar
    bipolar_data = []
    bipolar_names = []

    pairs = [
    ("CECOG_LF_1___01", "CECOG_LF_1___02"),
    ("CECOG_LF_1___02", "CECOG_LF_1___03"),
    ("CECOG_LF_1___03", "CECOG_LF_1___04"),
    ]

    for ch1, ch2 in pairs:

        x1 = xr_data.sel(chan=ch1).data
        x2 = xr_data.sel(chan=ch2).data

        bipolar_data.append(x1 - x2)
        bipolar_names.append(f"{ch1}-{ch2}")

    bipolar_data = np.array(bipolar_data)

    xr_bipolar = xr.DataArray(
        bipolar_data,
        dims=("chan", "time"),
        coords={
            "chan": bipolar_names,
            "time": xr_data["time"].values
        }
    )

    plot_xr_channels(xr_bipolar, chan_start=0, n_chan=3, target_srate=100, srate=500, offset=4, include_resp=True, norm_resp="zscore", vlines=None, title=None,)


    #### compare
    chan_i = 2

    plt.plot(xr_av_ref.data[chan_i], label='av_ref')
    plt.plot(xr_bipolar.data[chan_i], label='bipolar')
    plt.plot(xr_data.data[chan_i], label='raw')
    plt.legend()
    plt.show()
