
from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *
from A_config.n03_O_patient_info import *

import plotly
import plotly.graph_objects as go


debug = False








################################
######## INSPECT DATA ########
################################



#sujet = sujet_list[2]
def visualize_sig(sujet, step_group=20):

    #### get data
    print('IMPORT')
    data, resp, co2, chanlist, localist = get_data_raw(sujet)

    #### extract length
    chunk_time_list = np.cumsum([_chunk_data.shape[-1] for _chunk_data in data])

    #### generate linear signals
    print('CONSTRUCT LINEAR SIG')
    data_linear = np.concat(data, axis=1)
    resp_linear = np.concat(resp)
    co2_linear = np.concat(co2)

    #### downsample
    time_vec = np.arange(data_linear.shape[-1])
    data_down = scipy.signal.resample_poly(data_linear, up=1, down=10, axis=1, padtype='line')
    respi_down = scipy.signal.resample_poly(resp_linear, up=1, down=10, padtype='line')
    co2_down = scipy.signal.resample_poly(co2_linear, up=1, down=10, padtype='line')
    time_down = scipy.signal.resample_poly(time_vec, up=1, down=10, padtype='line')

    if debug:

        for chan_i in range(10):
            plt.plot(data_down[chan_i])
            plt.show()

        plt.plot(time_down)
        plt.show()

    #### plot
    print("PLOT")
    changroup_plot_list = np.arange(len(chanlist), step=step_group)

    #changroup_plot = changroup_plot_list[3]
    for changroup_plot in changroup_plot_list:

        fig = plotly.graph_objects.Figure()

            #### respi
        fig.add_trace(
                go.Scatter(
                    x=time_down,
                    y=scipy.stats.zscore(respi_down),
                    mode="lines",
                    name='respi'
                )
            )
        
            #### CO2
        y = scipy.stats.zscore(co2_down) - (1 * 4)
        
        fig.add_trace(
                go.Scatter(
                    x=time_down,
                    y=y,
                    mode="lines",
                    name='CO2'
                )
            )

            #### DATA
        chan_sel = np.arange(changroup_plot, changroup_plot+step_group)
        for chan_i_i, chan_i in enumerate(chan_sel):

            y = scipy.stats.zscore(data_down[chan_i]) - ((chan_i_i+2) * 4)

            fig.add_trace(
                go.Scatter(
                    x=time_down,
                    y=y,
                    mode="lines",
                    name=f"{localist[chan_i]} {chanlist[chan_i]}"
                )
            )

        for line in chunk_time_list:
            fig.add_vline(
                x=line,
                line_width=2,
                line_dash="dash",
                line_color="red"
            )

        fig.update_layout(
            template="simple_white",
            title=f"{sujet}",
            xaxis_title="Samples",
            yaxis_title="Amplitude (stacked)",
            height=800
        )

        fig.show()















################################
######## EXECUTE ########
################################

if __name__ == '__main__':

    sujet = 'NS217'
    sujet = 'LH018'

    visualize_sig(sujet, step_group=20)



