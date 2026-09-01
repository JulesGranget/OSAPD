function [epochs, timeVectors, timeTracker] = BU_EpochDataTrials(data, oldfs, newfs, trial_starts, trial_stops, varargin)    
% BU_EpochDataTrials    
% - If TDT_trial_lengths is provided (varargin{1}), treats input as SmartLab data    
%   and interpolates each SL trial to match TDT sample counts.    
% - If not provided, treats input as TDT data and returns sliced trials without interpolation.    
% - For eyetracker data: crops to TDT length if longer, pads with NaN if shorter if 'true' passed (NO EXTRAPOLATION)    
%    
% Inputs:    
%   data            - [chans x samples] or vector    
%   oldfs           - original sampling rate of 'data'    
%   newfs           - desired sampling rate after resampling    
%   trial_starts    - vector of trial start times (seconds)    
%   trial_stops     - vector of trial stop times (seconds)    
%   varargin{1}     - optional: TDT_trial_lengths (seconds)    
%   varargin{2}     - optional: isET (boolean, default=false) - if true, pad with NaN instead of extrapolating  
%    
% Outputs:    
%   epochs          - cell array of trials [chans x samples]    
%   timeVectors     - cell array of time vectors for each trial    
%   timeTracker     - [nTrials x 2] start/stop time in seconds

nTrials = length(trial_starts);    
epochs = cell(nTrials,1);    
timeVectors = cell(nTrials,1);    
timeTracker = zeros(nTrials,2);

is1D = isvector(data);    
if is1D    
    data = data(:);  % column [samples x 1]    
end

if oldfs ~= newfs    
    [p,q] = rat(newfs/oldfs);    
    if is1D    
        data_resampled = resample(data, p, q);    
    else    
        data_resampled = resample(data.', p, q).';    
    end    
else    
    data_resampled = data;    
end

% Sample dimension length    
if is1D    
    nSamples = numel(data_resampled);    
else     
    nSamples = size(data_resampled, 2);    
end

% Parse varargin  
doInterpToTDT = false;  
TDT_trial_lengths = [];  
isET = false;  % default: NOT eye tracker data

if ~isempty(varargin)    
    TDT_trial_lengths = varargin{1};  
    doInterpToTDT = true;  
      
    if numel(TDT_trial_lengths) ~= nTrials    
        error('TDT_trial_lengths must match number of trials.');    
    end  
      
    % Check for second varargin (isET flag)  
    if length(varargin) >= 2  
        isET = varargin{2};  
        if ~islogical(isET) && ~isnumeric(isET)  
            error('varargin{2} (isET) must be a boolean (true/false or 1/0)');  
        end  
        isET = logical(isET);  % ensure it's boolean  
    end  
end

for k = 1:nTrials

    start_idx = floor(trial_starts(k) * newfs) + 1;    
    start_idx = max(1, min(start_idx, nSamples));

    dur_sec = max(0, trial_stops(k) - trial_starts(k));    
    N = max(0, ceil(dur_sec * newfs));

    % Stop index, clamp    
    stop_idx = start_idx + N - 1;    
    stop_idx = min(stop_idx, nSamples);

    % Extract trial    
    if is1D    
        trial = data_resampled(start_idx:stop_idx);    
        SL_N_actual = length(trial);    
    else    
        trial = data_resampled(:, start_idx:stop_idx);    
        SL_N_actual = size(trial, 2);    
    end

    t_SL   = (0:SL_N_actual - 1) / newfs;    
        
    if doInterpToTDT    
        % Interpolate / crop trial to TDT trial length    
        TDT_dur = TDT_trial_lengths(k);    
        TDT_N   = max(0, ceil(TDT_dur * newfs));       
        t_TDT   = (0:TDT_N - 1) / newfs;    
        
        curr_N = length(t_SL);  % actual trial length    
        
        if curr_N < TDT_N    
            % Data is SHORTER than TDT  
                
            if isET    
                % EYE TRACKER: PAD with NaN instead of extrapolating    
                if is1D    
                    epochs{k} = [trial(:); NaN(TDT_N - curr_N, 1)];    
                else    
                    epochs{k} = [trial, NaN(size(trial,1), TDT_N - curr_N)];    
                end    
            else    

                    
                if is1D    
                    temp = interp1(t_SL, trial, t_TDT, 'linear', 'extrap');    
                    epochs{k} = temp(:);    
                else    
                    trial_interp = zeros(size(trial,1), numel(t_TDT));    
                    for ch = 1:size(trial,1)    
                        trial_interp(ch,:) = interp1(t_SL, trial(ch,:), t_TDT, 'linear', 'extrap');    
                    end    
                    epochs{k} = trial_interp;    
                end    
            end    
            
        elseif curr_N > TDT_N    
            % Data is LONGER than TDT - CROP to TDT length    
                
            if is1D    
                epochs{k} = trial(1:TDT_N);    
            else    
                epochs{k} = trial(:, 1:TDT_N);    
            end    
        
        else    
            % Same length, do nothing    
            epochs{k} = trial;    
        end    
        
        timeVectors{k}   = t_TDT;    
        timeTracker(k,:) = [trial_starts(k), trial_starts(k) + TDT_dur];    
    else    
        % TDT data path: no interpolation, just return the slice    
        epochs{k}        = trial;    
        timeVectors{k}   = t_SL;    
        timeTracker(k,:) = [trial_starts(k), trial_stops(k)];    
    end    
end    
end  
