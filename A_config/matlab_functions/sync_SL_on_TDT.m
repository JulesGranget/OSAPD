function data_sync = sync_SL_on_TDT(data_tosync, peaks_TDT, peaks_tosync, debug_token)

    if debug_token
        figure;
        plot(data_tosync(1,:));
    end

    data_sync = [];

    for data_i = 1:size(data_tosync, 1)

        data_sync_onechan = [];
        
        for seg_i = 1:(size(peaks_tosync,2)-1)
        
            n_tosync  = peaks_tosync(seg_i+1) - peaks_tosync(seg_i);
            n_TDT = peaks_TDT(seg_i+1) - peaks_TDT(seg_i);
        
            time_vec_tosync = linspace(0, 1, n_tosync);
            time_vec_TDT = linspace(0, 1, n_TDT);
                    
            if debug_token
                figure;
                plot(time_vec_SL); hold on;
                plot(time_vec_TDT);
                legend({'SL', 'TDT'});
            end
                    
            data = data_tosync(data_i, peaks_tosync(seg_i)+1:peaks_tosync(seg_i+1));
    
            data_synced = interp1( ...
                time_vec_tosync, ...
                data, ...
                time_vec_TDT, ...
                'linear', ...
                'extrap' ...
            );

            data_sync_onechan = [data_sync_onechan, data_synced];
        
                        
            if debug_token
                sig_plot = data_synced(1,:);
        
                figure;
                plot(1:numel(sig_plot), sig_plot); hold on;
                xline(peaks_TDT, 'r');
        
                figure;
                plot(1:numel(sig_plot), sig_plot); hold on;
                xline(peaks_tosync, 'r');
    
            end

        end

        data_sync = [data_sync; data_sync_onechan];

    end
