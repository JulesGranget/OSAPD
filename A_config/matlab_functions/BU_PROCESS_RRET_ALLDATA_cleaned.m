clear all; close all
%% JA

restoredefaultpath
addpath("/home/david/MATLAB Add-Ons/Collections/FieldTrip")
ft_defaults
addpath /home/david/Downloads/scripts
addpath /home/david/Downloads/scripts/fieldtrip-20230118/external/dmlt/external/murphy/KPMtools
addpath /home/david/Documents/RRET/BU/fct_jose/
addpath /home/david/Documents/EPIPE/Matlab/EXTERNAL/TDT/MatlabSDK/TDTSDK/TDTbin2mat
addpath(genpath('/home/david/Documents/RRET/BU/EPIPROC-master'))
addpath('/home/david/Documents/RRET/BU/FT_private/')
addpath('/home/david/Documents/RRET/BU')
addpath '/home/david/Documents/RRET/latestscripts_JA'
addpath('/home/david/Documents/RRET/samruddhi')

%% -----------------  SELECT PATIENT -----------------
pat = 'LH018'; %NS215, %NS213; %NS131_02; %'NS211';%'LH013';%'NS200';%'NS057_3';%'LH013';%'NS189';%'LH014_3';%'NS191';%'LH014_3'; 'NS057_3'
run 'LIST_filenames_RRET.m';

%RRET
%pat = 'NS202';
%run '/Volumes/F_BBio/RRET/LIST_filenames_RRET.m';

%% ----------------- PRE-PROCESSING STEPS ---------------------------------------------------------------------

% Steps 1-6: Pre-processing 
% Output structure contains:
% 1) tdt/natus data: eeg, resp, cranial, emg, plus anag chans (ttl,resp,...)
% 2) smartlab data : pressure, flow,
% Output is an 'epochs' structure with all synched data
% (all data synchronized and epoched with same fs)

%% ----------------------
Step1   = 0;
% Import TDT/Natus data and save as an ecog variable using specs from paramsfile

%% ----------------------
Step2   = 0;
% Visual inspection - Rewrite ecog fields

%% ----------------------
Step3   = 1;
% Filter, dowsample, save new vv ('..._epoch.m')

%% ----------------------
Step4 = 1;

%% ----------------- ANALYSES ------------------------------------------------------------------
% AFTER COMPUTING THIS STEP, JUMP TO ALLPATS_RRST_RESULTS.m

Find_ChanswithRRBOs = 0; 
% average all inhs and find chans with significant respiratory responses (ERPs)
% two steps: 1) epoching RRBOs 2) computing stats
% input data should be 'epochs' structure with TDT/SL fields & inspected inh onset times
% can include ALL inhs (loaded and non-loaded)
% epoch time-window: 2s before inh onset and 3s after
% include ALL trials (bad ones also as we need as many inhs as possible for RRBO analyses)

% 1) align epochs to inhonsets (not sink pulses) - this is done in P100 - NO NEED TO ALIGN IN RRST/RRET
% 3) add fields to 'epochs' - options are filter, bipolarize, discard inhns with IEDs.
% 4) add field 'epochs.RBBO': info about chans with *ERPs, p-vals, ...
% 5) add subfields to 'epochs.RBBO': spectrum, itc, p-vals


%% --------------------------------------------------------------------------------------------

if Step1 %% use params.m; discard bad chans; assign trigger chans;
    % run preproc
    % go to subject folder (projects/resp/respcceps/ns141/N141_preproc.m and
    if strcmp(params.recsystem,'TDT')
        [ecog] = TDT2ecog(params.paramsfile);
    elseif strcmp(params.recsystem,'NATUS')
        %[ecog] = BU_Natus2mat_quick(params);
        [ecog] = NATUS2ecog(params);
    end
end

%% --------------------------------------------------------------------------------------------
if Step2
    %% inspects signal - visualize iEEG chan
    %% creates _epoch vv but does not apply any filtersStep

    %% Or load the ecog variable if step 1) was done already
    load([params.directoryOUT params.filename '_ecog.mat'])
    ecog.params = params;%%replace params in case it has not been updated

    %% 2) Check for bad channels and update xls
    % Look at EEG, write bad channels to Xls-sheet
    cfg =[];
    cfg.channel = 'all';%'LIa*';%{'LFo*'};%'all';
    %ttl = ecog.ttl/max(ecog.ttl);voltagethreshold = 0.2;
    %[timestamps] = BU_find_timestampsTTLs(ttl,voltagethreshold,ecog.ftrip.fsample)
    %timestamps(1)=[];

    eeg = ecog.ftrip.trial;%ecog.ftrip.trial{1};
    fs  = ecog.ftrip.fsample;

    %ecog_databrowser(ecog,cfg)

    %% if eeg not looking good try
    %% remove notch
    selchan = 20;%just for visualization
    if iscell(eeg); eeg = eeg{1}; end % JA added

    eeg = fillmissing(eeg, 'linear', 2);   %sam added to handle NaNs in eeg
    eeg = double(eeg); % sam added due to some precision errors during filtering (inf/nans)

    [eeg1] = BU_filternoise_harmonics(eeg,fs,[],selchan);
    [eeg1] = BU_filternoise_harmonics(eeg,fs,params.filter_notch,selchan);

    %% filter eeg (bandpas)
    n = 2; Wn = [ecog.params.filter_bandpass]/(fs/2);[b,a] = butter(n,Wn);%jos
    eeg2 = filtfilt(b,a,double(eeg1'))';

    figure
    plot(demean(eeg(5,fs:fs*10)),'b');hold on
    plot(demean(eeg1(5,fs:fs*10)),'m');hold on
    plot(demean(eeg2(5,fs:fs*10)),'r');hold on

    %% collect (temp)
    if ~iscell(ecog.ftrip.trial) % JA
        ecog.ftrip.trial = {eeg2}; 
    end
    if ~iscell(ecog.ftrip.time) % JA
        ecog.ftrip.time = {ecog.ftrip.time}; 
    end
    
    ecog.ftrip.trial{1} = eeg2; %% only for visualization (this get
    % re-written bellow) % commented by JA
    ecog.params.filter_notch = params.filter_notch;%% = 1 % annotate in params file what you wanna do in step3 (notch filtering)
    ecog.params.filter_bandpass = params.filter_bandpass; %% = [0.5 200];% annotate in params file what you wanna do in step3 (bandpass filtering)
    %%visualise
    ecog_databrowser(ecog,cfg)

    %%save timestamps using the browser to get start and end of sz (save to
    %%params file also and ecog vv.


    %% save raw signal in _ecog file %% (filters should be applied in next step)
    ecog.ftrip.trial{1} = eeg;

    %% Look at REF (skull) electrodes and write selected ones in params file (params.goodchancran=[2 3])
    if isfield(ecog.params, 'cranial') && ~isempty (ecog.params.cranial)
        close all
        for chan = 1:size(ecog.cranial.trial,1)
            figure;plot(ecog.cranial.trial(chan,:))
        end
        disp('write to params file good ref/skull elec');pause  %%if no goog >>> params.goodchancran=[];
        if isempty(params.goodchancran)
            ecog.cranial=[];
        end
    end

    % If you added bad (or SOZ, spikey, out) chans to xls. Read xls in again and save ecog variable.
    % (this will overwrite several fields in the ecog structure)
    [labelfile,ecog] = BU_read_labelfile(params.labelfile,ecog);
    save([params.directoryOUT params.filename '_ecog.mat'],'ecog');

    % you can also use [labelfile,ecog] = read_labelfile(params.labelfile,ecog,1); if you
    % don't want to be asked if you want to overwrite.

    %% 4) Find bad channels using PSD (not mandatory)
    find_bad_chans(ecog);
    save([params.directoryOUT params.filename '_ecog.mat'],'ecog','-v7.3');

    %% if you want to look at PSD (FFT) at any point
    cfg = [];
    cfg.plot = 1;
    cfg.T = 3; % epoch size in seconds
    ecog  = ecog_psd(ecog,cfg);
    %save([params.directoryOUT params.filename '_ecog.mat'],'ecog');

    cfg=[];
    cfg.chans = 'all';
    plot_ecog_psd_ft(ecog)
end

%% ------------------------ STEP 3 ------------------------------------------------
%% discard bad electrodes
%% filter (band-pass,notch) and downsamples
%% make all signals equal length (iEEG, resp, ekg,...)
%% adds smartlab data
%% find timestamps (done in next step)
if Step3

    clear Step*

    load([params.directoryOUT params.filename '_ecog.mat'])
    ecog.params = params;%%replace params in case it has not been updated
    %%%save([params.directoryOUT params.filename '_ecog.mat'],'ecog');%rewrite if neccessry

    plotting = 1;

    %% discard electrode contacts: out/bad and WM (based on excell sheet)
    % idxBAD = find(~isnan(ecog.labelimport.bad));
    % idxOUT = find(~isnan(ecog.labelimport.out));
    % %idxWM  = find(ecog.labelimport.labels_FSurf_GM~=1);%%exclude WM & CSF electrodes
    % idxWM  = [];%find(ecog.labelimport.labels_FSurf_GM==2);%%exclude CSF electrodes
    % temp = [];
    % if ~isempty(idxBAD);temp = [temp; idxBAD];end
    % if ~isempty(idxOUT);temp = [temp; idxOUT];end
    % if ~isempty(idxWM); temp = [temp; idxWM];end
    % 
    % idxOFF = unique(temp);
    % idxOFFlab = ecog.labelimport.labels(idxOFF);
    % idxOFFlab = strrep(idxOFFlab,' ','');%%remove epmty
    % 
    % idxOFFf = [];
    % for n = 1:length(ecog.ftrip.label)
    %     x = strcmp(ecog.ftrip.label(n),idxOFFlab);
    %     if sum(x)==1
    %         idxOFFf = [idxOFFf n];
    %     end
    % end
    % 
    % idxIN        = ~ismember(1:length(ecog.ftrip.label),idxOFFf);
    % % JA added cell conversion below
    % if ~iscell(ecog.ftrip.trial)
    %     ecog.ftrip.trial = {ecog.ftrip.trial};
    % end
    % if ~iscell(ecog.ftrip.time)
    %     ecog.ftrip.time = {ecog.ftrip.time};
    % end
    % 
    % eeg          = ecog.ftrip.trial{1}(idxIN,:);
    % eeg_raw      = eeg; % JA: to save unfiltered EEG later
    % labels       = ecog.ftrip.label(idxIN);
    % labels_FSurf = ecog.ftrip.label_FSurf(idxIN);
    % clear idxIN idxBAD idxOUT idxOFF idxWM idxOFFf n x
    % fs = ecog.ftrip.fsample;

    %% ------ reject bad skull elecs ------
    % if isfield (params, 'goodchancran') && ~isempty (params.goodchancran)
    % 
    %     if iscell (ecog.cranial.trial)
    %         cranial = ecog.cranial.trial{1}(params.goodchancran,:);
    %     else
    %         cranial = ecog.cranial.trial(params.goodchancran,:);
    %     end
    % 
    %     cranial_label = ecog.cranial.label(params.goodchancran);
    % end

    %% ------ EEG data ------

    % -- band-pass filtering
    eeg          = ecog.ftrip.trial{1};
    eeg_raw      = eeg; % JA: to save unfiltered EEG later
    labels       = ecog.ftrip.label;
    labels_FSurf = ecog.ftrip.label_FSurf;
    clear idxIN idxBAD idxOUT idxOFF idxWM idxOFFf n x
    fs = ecog.ftrip.fsample;

    if ~isempty (params.filter_bandpass) %%filter eeg (bandpass)

        figure;plot(eeg(1,fs:fs*10));hold on;
        n = 2; Wn = [ecog.params.filter_bandpass]/(fs/2);[b,a] = butter(n,Wn);%jos
        eegF1 = filtfilt(b,a,double(eeg)')';
        plot(eegF1(1,fs:fs*10),'r');hold on;
        filter_band = ecog.params.filter_bandpass;
        if filter_band(1) == 0;n = 2; wn = filter_band(2)/(fs/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
        if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(fs/2);[b,a] = butter(n,wn);end%%%band-pass

        % -- use FiltFiltM
        current_directory = pwd;
        % temp_directory = '/Users/jherreroru/Documents/gitmatlab/fct';
        % cd(temp_directory);
        eegF2 = FiltFiltM(b,a,double(eeg),2);
        % cd(current_directory);
        plot(eegF2(1,fs:fs*10),'g');hold on;
        eeg = eegF2;% **
        clear eegF1 eegF2
    end

    if ~isempty(params.filter_notch)
        selchan = 1;%chan to plot -filter applies to all chans
        filternoise_harmonics = [params.filter_notch];
        [eeg] = BU_filternoise_harmonics(eeg,fs,filternoise_harmonics,selchan);
    end

    % --- downsampling eeg ---
    downsampling = 1;
    if downsampling
        newfs = 500;oldfs = fs;
        [eeg,fs] = BU_ResampleData(eeg',newfs,oldfs,params.recsystem);
        % %JA made own version, should replace ***
        clear newfs oldfs
    end

    % --- demean eeg
    %%%eeg = single(demean(eeg)); %not neccesary

    %% ------ EMG data ------
    if isfield(ecog,'emg') && ~isempty(ecog.emg)
        emg = ecog.emg.trial{:}; %plot(emg(1,1:20000))
        emg = double(fillmissing(emg, 'linear', 2)); % sam added to take care of NaNs

        if ~isempty (params.filter_bandpass)
            % band-pass (butter)
            filter_band = params.filter_bandpass;%[0 200];
            if filter_band(1) == 0;n = 2; wn = filter_band(2)/(fs/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
            if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(fs/2);[b,a] = butter(n,wn);end%%%band-pass

            current_directory = pwd;% use FiltFiltM
            temp_directory = '/Users/jherreroru/Documents/gitmatlab/fct';
            % cd(temp_directory);
            emgF = FiltFiltM(b,a,emg,2);
            % cd(current_directory);
            plotting = 1;
            if plotting;plot(emg(1,1:2000));hold on; plot(emgF(1,1:2000),'r');end
            emg = emgF;clear emgF;
        end

        % notch filtering
        if ~isempty(params.filter_notch)
            selchan = 1;%chan to plot -filter applies to all chans
            filternoise_harmonics = [params.filter_notch];
            
            if isfield (ecog.emg, 'fs')
                emg_fs = ecog.emg.fs;
            end

            if isfield (ecog.emg, 'fsample')
                emg_fs = ecog.emg.fsample;
            end
            
            [emg] = BU_filternoise_harmonics(emg,emg_fs,filternoise_harmonics,selchan);
        end

        % downsampling emg
        if emg_fs ~= fs && downsampling
            newfs = fs;oldfs = emg_fs;
            [emg] = BU_ResampleData(double(emg)',newfs,oldfs,params.recsystem);
            clear newfs oldfs
        end

        % check lengths
        if length(emg) > length(eeg)
            emg = emg(:,1:legth(eeg));
            disp('check emg length');pause
        end

        % current_directory = pwd;
        % temp_directory = '/Users/jherreroru/Documents/GitHub/EPIPE/Matlab/utilities';
        % cd(temp_directory);
        if size(emg,1) > 1
            emg = single(demean(emg,2));
        else
            emg = single(demean(emg,1));
        end
        emg_label = ecog.emg.label; %plot(emg(1,1:20000))
        % cd(current_directory);
    end

    %% ---- SCALP data - filtering & downsampling -------
    if isfield(ecog,'scalp') && ~isempty(ecog.scalp)
        scalp = ecog.scalp.trial{:}; %plot(scalp(1,1:20000))

        % band-pass filter
        if ~isempty(params.filter_bandpass)
            if filter_band(1) == 0;n = 2; wn = filter_band(2)/(fs/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
            if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(fs/2);[b,a] = butter(n,wn);end%%%band-pass
            scalpF = FiltFiltM(b,a,scalp,2);
            plotting = 1;
            if plotting;plot(scalp(1,1:2000));hold on; plot(scalpF(1,1:2000),'r');end
            scalp = scalpF;clear scalpF;
        end

        % notch removal
        if ~isempty(params.filter_notch)
            selchan = 1;%chan to plot -filter applies to all chans
            filternoise_harmonics = [params.filter_notch];
            [scalp] = BU_filternoise_harmonics(scalp,ecog.scalp.fs,filternoise_harmonics,selchan);
        end

        % downsampling
        if ecog.scalp.fsample ~= fs && downsampling
            newfs = fs; oldfs = ecog.scalp.fsample;
            [scalp] = BU_ResampleData(scalp',newfs,oldfs,params.recsystem);
            clear newfs oldfs
        end
        if length(scalp) > length(eeg)
            scalp = scalp(:,1:legth(eeg));
        end
        scalp = single(demean(scalp,2));
        scalp_label = ecog.scalp.label; %plot(emg(1,1:20000))
    end

    %% ------ CRANIAL data ------
    if isfield(ecog,'cranial') && isfield(ecog,'params') && ~isempty(ecog.cranial) && ~isempty(params.goodchancran)
        if isfield(ecog.cranial, 'fsample')
            fsnow = ecog.cranial.fsample;
        elseif isfield(ecog.cranial, 'fs');
            fsnow = ecog.cranial.fs;
        else
            fsnow = fs;
        end

        % band-pass
        if ~isempty(params.filter_bandpass)
            filter_band = params.filter_bandpass;
            if filter_band(1) == 0;n = 2; wn = filter_band(2)/(fsnow/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
            if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(fsnow/2);[b,a] = butter(n,wn);end%%%band-pass
            %cranialF = FiltFiltM(b,a,cranial,2);
            cranialF = filtfilt(b,a,double(cranial)')';
            plotting = 1; if plotting;close all;plot(cranial(1,1:4000));hold on; plot(cranialF(1,1:4000),'r');end
            cranial = demean(cranialF,2);
            clear cranialF;
        end

        % notch removal
        if ~isempty(params.filter_notch)
            selchan = 1;%chan to plot
            
            if iscell(ecog.cranial.trial)
                [cranial] = BU_filternoise_harmonics(ecog.cranial.trial{1},fsnow,params.filter_notch,selchan);
            else
                [cranial] = BU_filternoise_harmonics(ecog.cranial.trial,fsnow,params.filter_notch,selchan);
            end

        end

        % downsampling cranial
        if fsnow ~= fs || downsamplingS == 1
            if fsnow ~= fs && downsampling
                newfs = fs;oldfs = fsnow;
                [cranial] = BU_ResampleData(double(cranial)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end
        end

        % check lengths
        if length(cranial)> length(eeg)
            cranial=cranial(:,1:length(eeg));
            disp('check length of cranial data');pause
        end
        cranial = single(demean(cranial,2));
        cranial_label; %plot(emg(1,1:20000))
    end

    %% -------- filtering & downsampling resp -------
    if isfield(ecog.params,'resp') && ~isempty(ecog.params.resp) && ~isempty(ecog.params.resp.chan)
        if isfield(ecog.resp, 'trial')
            resp = ecog.resp.trial{:};
            resp_fs = ecog.resp.fs;
        else
            resp = ecog.resp;
            disp("is resp_fs = eeg_fs?");pause
            resp_fs = ecog.ftrip.fsample;
        end %plot(resp(1,1:20000))

        filter_band=[0 2];
        if filter_band(1) == 0;n = 2; wn = filter_band(2)/(resp_fs/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
        if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(resp_fs/2);[b,a] = butter(n,wn);end%%%band-pass


        current_directory = pwd;% use FiltFiltM
        % temp_directory = '/Users/jherreroru/Documents/gitmatlab/fct';
        % cd(temp_directory);
        respF = FiltFiltM(b,a,resp,2);
        % cd(current_directory);
        plotting=1;
        if plotting;plot(resp(1,1:end));hold on; plot(respF(1,1:end),'r','Linewidth',4);end
        resp = respF;clear respF;close all hidden

        % --- downsampling
        if ecog.resp.fs ~= fs && downsampling
            newfs = fs; oldfs = ecog.resp.fs;
            [resp] = BU_ResampleData(resp',newfs,oldfs,params.recsystem); % JA changed resp to resp'
            clear newfs oldfs
        end

        if length(resp) > length(eeg)
            disp("check lenghts");pause
            resp = resp(:,1:length(eeg));
        end

        resp = single(demean(resp,2));
    end

    %% ---------- take resp from ana -------
    if isfield(ecog,'analog') && isfield(ecog.analog,'label')
        if any(strcmp(ecog.analog.label,'resp'))
            idx = find(strcmp(ecog.analog.label,'resp'));
            % idx = find(cellfun(@isempty,x)==0);
            if ~isempty(idx)
                resp_ana = ecog.analog.trial{1}(idx,:);
                resp_ana_fs = ecog.analog.fs;
            end
            if resp_ana_fs ~= fs
                filter_band=[0 2];
                if filter_band(1) == 0;n = 2; wn = filter_band(2)/(resp_ana_fs/2);[b,a] = butter(n,wn,'low');end%%%low-pass filter <200Hz
                if filter_band(1) ~= 0;n = 2; wn = [filter_band]/(resp_ana_fs/2);[b,a] = butter(n,wn);end%%%band-pass
                resp_anaF = FiltFiltM(b,a,resp_ana,2);
                plotting=1;
                %if plotting;plot(resp_ana(1,1:200000));hold on; plot(resp_anaF(1,1:200000),'r');end
                if plotting;plot(resp_ana(1,1:183808));hold on; plot(resp_anaF(1,1:183808),'r');end
                resp_ana=resp_anaF;clear resp_anaF;
                %%downsampling
                downsampling = 1;
                if downsampling
                    dataf = 500;%new fs
                    fs_data = resp_ana_fs;%old fs
                    % [upsamp_data,downsamp_data] = rat(dataf/fs_data,1e-10);
                    % resp_ana = resample(double(resp_aresp_anana)', upsamp_data, downsamp_data)';
                    resp_ana = BU_ResampleData(resp_ana', dataf, fs_data, params.recsystem);
                    clear dataf upsamp_data downsamp_data fs_data wn n a b
                end
            end
            if length(resp_ana)> length(eeg)
                disp("check lengths");pause
                resp_ana = resp_ana(:,1:length(eeg));
            end
            resp_ana = single(demean(resp_ana));
            resp = resp_ana; %samruddhi added so resp is taken from analog channel
            %resp = resp(2,:); %plot(emg(1,1:20000))
        end
        ecog.resp = resp; %samruddhi added; to ensure that analog resp channel data is saved in ecog. No ecog.resp created otherwise
    end

    %% get only 1 chan for resp (inh going downward
    % if exist('resp','var') && ~isempty(resp) && size(resp,1) > 1
    %     % figure;
    %     % plot(resp(1,:),'r');hold on
    %     % plot(resp(2,:),'k')
    % 
    %     % Check the user response
    %     user_response = input('Select ONLY ONE resp channel: 1 or 2: ', 's');
    % 
    %     if strcmpi(user_response, '1')
    %         resp = resp(1,:);
    %     elseif strcmpi(user_response, '2')
    %         resp = resp(2,:);
    %     end
    %     close all hidden; clear user_response
    % end
    %% sam added for gsm mixer - output CO2 flow
    if isfield(ecog, 'analog') && isfield(ecog.analog, 'label')
        if length(ecog.analog.label) >= 3
            gsm_mixer_label = ecog.analog.label{3}; % Label from channel 3
            
            % Check if the label matches the intended use case ('gsm mixer', adjust as needed)
            if any(strcmpi(ecog.analog.label, gsm_mixer_label))
                % Extract the data from the specified channel
                gsm_mixer_ana = ecog.analog.trial{1}(3, :);
                gsm_mixer_fs = ecog.analog.fs;
    
                % Resample if necessary
                if gsm_mixer_fs ~= fs
                    filter_band = [0 2];

                    % Filtering
                    n = 2;
                    if filter_band(1) == 0
                        wn = filter_band(2) / (gsm_mixer_fs / 2);
                        [b, a] = butter(n, wn, 'low');
                    else
                        wn = filter_band / (gsm_mixer_fs / 2);
                        [b, a] = butter(n, wn);
                    end

                    gsm_mixer_anaF = FiltFiltM(b, a, gsm_mixer_ana, 2);

                    gsm_mixer_ana = gsm_mixer_anaF;
                    clear gsm_mixer_anaF;
    
                    % Downsampling
                    downsampling = 1;
                    if downsampling
                        dataf = 500; % New fs
                        fs_data = gsm_mixer_fs; % Old fs
                        % [upsamp_data, downsamp_data] = rat(dataf / fs_data, 1e-10);
                        % gsm_mixer_ana = resample(double(gsm_mixer_ana)', upsamp_data, downsamp_data)';
                        gsm_mixer_ana = BU_ResampleData(gsm_mixer_ana', dataf, fs_data, params.recsystem);

                        clear dataf upsamp_data downsamp_data fs_data wn n a b
                    end
                end
    
                % Ensure the length matches the EEG data or other baseline data
                if length(gsm_mixer_ana) > length(eeg)
                    disp("check lengths");
                    pause;
                    gsm_mixer_ana = gsm_mixer_ana(:, 1:length(eeg));
                end
    
                gsm_mixer_ana = single(demean(gsm_mixer_ana));
    
                ecog.gsm_mixer = gsm_mixer_ana;
            end
        end
    end

    %% ---------------- SMARTLAB --> TDT  ------------------------
    % connections can be:
    % 1) direct --> from SL digital inputs to TDT analogs (SINK pulses - 0ms delay)
    % 2) direct --> from SL digital outputs to TDT analogs (~40ms delay)
    % 3) indirect --> via usb  Instacall (~2s delay due to Instacall/Windows delays)
    %    from SL to TDT-analogs (via Instacal usb)
    %    correct time delay in TDT analog chans 
    % NOTE: this time delays DO NOT apply to the SL XML.CSV files

    %only commeted for NS221 with no TTL
    %for VCCL BASELINE

    if isfield(params.smartlab, 'tdt') && ~isempty(params.smartlab.tdt)

        close all hidden

        SmartLab = [];
        Ssignal  = [];

        correctdelay = 1;

        %% Sinking pulses -- Important!!
        % tdt channel used for sinking pulses (usually ana#1 diginput -direct)

        if isfield(params.smartlab,'sinkchan') && ~isempty(params.smartlab.sinkchan)

            sinkpulses = [];

            Schan = params.smartlab.sinkchan(1);
            Stype = params.analog.label{Schan};

            Ssignal = demean(ecog.analog.trial{1}(Schan,:));% 'demean' brings signal to zero
            Sfs = ecog.analog.fs;
            Slabel = ecog.analog.label{Schan};

            % check signal amplitude and rescale if needed
            % in natus is usually very small

            %if max(abs(Ssignal)) < 1
            Ssignal = Ssignal/max(abs(Ssignal));%range from 0-1
            %end

            figure;subplot(211);
            plot(Ssignal);hold on;

            data = [ 0 diff(Ssignal)];%this tends to normalize but not always (use below code if no good)
            plot(data,'r');hold on

            disp('if ttls are inverted then uncomment below this')
            pause
            % 
            % data = Ssignal*-1;
            % data = data-median(data);
            % clf;plot(data,'r');hold on

            % digitize signal: find timestamps -
            % set some params (threshold & min distance)

            thlim = 0.1; %2 % Set threshold
            mpdlim = round(Sfs);% Set your minimum peak distance here

            data_diff = [0, diff(data >=  thlim)];  % Take the difference of the thresholded data

            pulse_starts = find(data_diff == 1);  % Indices where the data first crosses the threshold

            LOCS = pulse_starts([true, diff(pulse_starts) >= mpdlim]);

            % KEEP COMMENTED
            % data_diff_fall = [0, diff(data >= thlim)];      
            % pulse_ends = find(data_diff_fall == -1);         % indices where it crosses from 1 -> 0
            % pulse_ends = pulse_ends([true, diff(pulse_ends) >= mpdlim]);  

            % visualize ends
            % xline(pulse_ends, 'b:'); % add blue lines for pulse ends

            % get pulse end times (seconds)
            % sinkpulses_ends = timenow(pulse_ends)';  % new variable with end times
            % -----------------------------------------------------------------------

            %%[PKS,LOCS] = findpeaks(data,'MinPeakHeight',thlim,'MinPeakDistance', mpdlim);
            %%[PKS,LOCS] = findpeaks(Ssignal,'MinPeakHeight',thlim,'MinPeakDistance', mpdlim);

            %KEEP COMMENTED

            % discard unwanted timestamps (based on params)
            if ~isempty(params.sampleON) && ~isempty(params.sampleOFF)
                idxnow = find (LOCS>params.sampleON & LOCS<params.sampleOFF);
                LOCS = LOCS(idxnow);
                %PKS  = PKS(idxnow);
            end
            % vline(LOCS,'r:');hline(thlim) % JA replaced with below
            xline(LOCS, 'r:');yline(thlim);
            title ('TDT sink pulses');

            % check timestamps not too close to each other
            subplot(212);
            plot(diff(LOCS),'r.');
            temp = round(min(diff(LOCS))/Sfs);
            title (['mindistance =' num2str(temp) 's']);

            disp(' Are SINK ONSETS accurate ?');pause

            timenow = (1:length(data))*1/Sfs;

            sinkpulses = timenow(LOCS)'; % get the times

            close all hidden

            %% if ttls are choppy, uncomment below
            %idx = [1 find(diff(LOCS)>Sfs)+1];%refractory period
            %vline(LOCS(idx),'k');
            %sinkpulses = timenow(LOCS(idx));

            % collect
            ecog.smartlab.tdt.sinkpulses = sinkpulses;

            clear temp timenow LOCS PKS
        end

        %%  ---- Inhalation onsets in TDTana#2 (sent from SL digout) ----

        if isfield (params.smartlab.tdt,'inhON') && ~isempty(params.smartlab.tdt.inhON)

            anachan = find(contains(params.analog.label,'inh(digout)')==1);%Smartlab digout#1 should be hooked to TDTana#2
            if ~isempty(anachan)
                Sdelay = params.smartlab.tdt.inhON(2);
            else
                anachan = find(contains(params.analog.label,'pressure')==1);
                Sdelay = 0.2;
            end

            Stype = params.analog.label{anachan};
            Schan = params.smartlab.tdt.inhON(1);
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs = ecog.analog.fs;
            Slabel = ecog.analog.label{Schan};

            % find Inhalation Onsets (in secs)
            inhonsets  = [];
            inhoffsets = [];

            figure; subplot(211)
            plot(Ssignal);hold on
            title('Original Data');

            % Ask user whether to improve snr or not
            choice = questdlg('Do you want to improve SNR ?', ...
                'clean signal ', ...
                'Yes', 'No', 'No');

            % Handle the response
            switch choice
                case 'Yes'
                    if ~exist('ft_preproc_bandpassfilter', 'file')
                        error('FieldTrip is not in the MATLAB path. Please add it using addpath.');
                    end

                    % Apply a band-pass filter using FieldTrip's ft_preproc_bandpassfilter function
                    low_freq = 10;
                    high_freq = 100;
                    filtered_data = ft_preproc_bandpassfilter(Ssignal, Sfs, [low_freq, high_freq]);

                    % Compute the first derivative
                    data_derivative = diff(filtered_data);
                    subplot(212)
                    plot(data_derivative);hold on
                    title('Cleaner data');

                    Ssignal = data_derivative./max(abs(data_derivative)); % collect normilized signal 0-1
                    clear data_derivative filtered_data low_freq high_freq

                case 'No'
                    disp('Signal not cleaned');
            end

            if contains(Stype,'ana') % anag signal (eg. pressure)
                inhUP = 0;minpeaktime = 0.75;minpeakamp = 0.1;
                [inhonsets,Ssignal] = BU_InhOnsets(Ssignal,Sfs,inhUP,minpeaktime,minpeakamp,Slabel);

            elseif contains(Stype,'dig') %ttl square waves

                x = [0 diff(Ssignal)];
                plot(x,'r')
                [PKS,LOCS] = findpeaks(x,'MinPeakHeight',.01,'MinPeakDistance',Sfs);

                % discard unwanted timestamps (based on params)
                if ~isempty(params.sampleON) && ~isempty(params.sampleOFF)
                    idxnow = find (LOCS>params.sampleON & LOCS<params.sampleOFF);
                    LOCS = LOCS(idxnow);
                    PKS  = PKS(idxnow);
                end
                vline(LOCS,'r:')
                inhonsets  = LOCS/Sfs; % in secs
                % correct SL-delay on time stamps
                if Sdelay ~= 0 && correctdelay
                    inhonsets = inhonsets - params.smartlab.tdt.inhON(2);%
                end

                % Ask user whether to get inh offsets
                choice = questdlg('Do you want to get inhalation offsets ?', ...
                    'inhalation offsets ', ...
                    'Yes', 'No', 'No');
                switch choice
                    case 'Yes'
                        [PKS,LOCS] = findpeaks(x*-1,'MinPeakHeight',.2);

                        % discard unwanted timestamps (based on params)
                        if ~isempty(params.sampleON) && ~isempty(params.sampleOFF)
                            idxnow = find (LOCS>params.sampleON & LOCS<params.sampleOFF);
                            LOCS = LOCS(idxnow);
                            PKS  = PKS(idxnow);
                        end

                        inhoffsets  = LOCS/Sfs; % in secs

                        %% correct for SL delays
                        if Sdelay ~= 0 && correctdelay
                            inhoffsets = inhoffsets - params.smartlab.tdt.inhON(2);%
                        end

                    case 'No'
                        disp('Not possible to calculate inh offsets with this signal');
                end

            end

            %% collect for SL delay on the continous signal

            if Sdelay ~= 0 && correctdelay

                % correct for SmartLab time delays
                delaysamples = round(params.smartlab.tdt.inhON(2)*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];

                % check matching lengths
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end

                Ssignal = Ssignal2; clear Ssignal2
            end

            % downsampling
            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            % collect
            ecog.smartlab.tdt.inhON = single(inhonsets);
            ecog.smartlab.tdt.inhOFF = single(inhoffsets);
            ecog.smartlab.tdt.fs = Sfs;

            if contains(Stype,'ana')
                ecog.smartlab.tdt.pressure = single(Ssignal);% collect - needs scaling
            elseif contains(Stype,'dig')
                ecog.smartlab.tdt.inhalation = single(Ssignal);% collect - needs scaling
            end

        end

        %% ---- balonONOFF (TDTana#)
        % as recorded in TDTana# (sent from SL digout) ----

        ecog.smartlab.tdt.balonONOFF = []; % digital output (ana2)
        if isfield (params.smartlab.tdt,'balonONOFF') && ~isempty (params.smartlab.tdt.balonONOFF)
            if params.smartlab.tdt.balonONOFF(2) == 0 && ~isnan(params.smartlab.tdt.balonONOFF(1))%real time (digital signal)
                Sdelay  = 0;
                Stype   = 'digital_output';
                Schan   = params.smartlab.tdt.balonONOFF(1);% chan used to determine balloon openings
                Ssignal = ecog.analog.trial{1}(Schan,:);
                Sfs     = ecog.analog.fs;
                Slabel  = ecog.analog.label{Schan};
            end
            Ssignal = Ssignal*-1;% invert it
            Ssignal = Ssignal - median(Ssignal);% zero it
            Ssignal = Ssignal./max(max((Ssignal))); % scale
            Ssignal = round(Ssignal);% round it - zeros show balloonOFF - ones show ballonON)

            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.balonONOFF = single(Ssignal);% collect
        end

        %% ---- Breathswitch TDTanag#
        % as recorded in TDTana# sent from SL digital input to TDTana#1 - direct NO delay) -----

        ecog.smartlab.tdt.breathswitch = []; % digital input
        if isfield (params.smartlab.tdt,'breathswitch') && ~isempty(params.smartlab.tdt.breathswitch)
            if params.smartlab.tdt.breathswitch(2) == 0 %real time (digital input)
                Sdelay  = 0;
                Stype   = 'digital_input';
                Schan   = params.smartlab.tdt.breathswitch(1);% chan with balloon openings
                Ssignal = ecog.analog.trial{1}(Schan,:);
                Sfs     = ecog.analog.fs;
                Slabel  = ecog.analog.label{Schan};
            end
            Ssignal = Ssignal*-1;% invert signal
            Ssignal = Ssignal - median(Ssignal);
            Ssignal = Ssignal./max(max((Ssignal)));
            Ssignal = round(Ssignal);% round it - zeros show balloonOFF - ones show ballonON)

            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.breathswitch = single(Ssignal);% collect
        end

        %% ---  Pressure (TDTanag#) -----
        % pressure as recorded in TDTanag# via usb
        % SL-->Instacall-->TDTana# [.2s delay]

        if isfield (params.smartlab.tdt,'pressure')  && ~isempty(params.smartlab.tdt.pressure) && ~isfield(ecog.smartlab.tdt,'pressure')

            % get delay params
            Sdelay  = params.smartlab.tdt.pressure(2);
            Stype   = 'analog';
            Schan   = params.smartlab.tdt.pressure(1);% analog chan were pressure was recorded
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs     = ecog.analog.fs;
            Slabel  = ecog.analog.label{Schan};

            % shift data backwards to correct for delay
            correctdelay = 1;

            if Sdelay ~= 0 && correctdelay
                delaysamples = round(Sdelay*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
                % check size of vector
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end
                Ssignal = Ssignal2; clear Ssignal2
            end

            % downsample signal
            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.pressure = single(Ssignal); % collect
        end

         %% ---  Flow (TDTanag#) -----
        % pressure as recorded in TDTanag# via usb
        % SL-->Instacall-->TDTana# [.2s delay]

        if isfield (params.smartlab.tdt,'flow')  && ~isempty(params.smartlab.tdt.flow) && ~isfield(ecog.smartlab.tdt,'flow')

            % get delay params
            Sdelay  = params.smartlab.tdt.flow(2);
            Stype   = 'analog';
            Schan   = params.smartlab.tdt.flow(1);% analog chan were pressure was recorded
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs     = ecog.analog.fs;
            Slabel  = ecog.analog.label{Schan};

            % shift data backwards to correct for delay
            if Sdelay ~= 0 && correctdelay
                delaysamples = round(Sdelay*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
                % check size of vector
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end
                Ssignal = Ssignal2; clear Ssignal2
            end

            % downsample signal
            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.flow = single(Ssignal); % collect
        end


        %% Volume
        % Volume as recorded on TDTanag# (SL-->Instacall-->TDTanag#)

        if isfield (params.smartlab.tdt,'volume')  && ~isempty(params.smartlab.tdt.volume) && ~isfield(ecog.smartlab.tdt,'volume')

            Sdelay  = params.smartlab.tdt.volume(2);
            Stype   = 'analog';
            Schan   = params.smartlab.tdt.volume(1);% analog chan were volume was recorded
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs     = ecog.analog.fs;
            Slabel  = ecog.analog.label{Schan};

            if Sdelay ~= 0 && correctdelay
                delaysamples = round(Sdelay*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
                % check size of vector
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end
                Ssignal = Ssignal2; clear Ssignal2
            end
            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.volume = single(Ssignal); % collect
        end

        %% ---  CO2 -----
        if isfield (params.smartlab.tdt,'co2')  && ~isempty(params.smartlab.tdt.co2)

            Sdelay  = params.smartlab.tdt.co2(2);
            Stype   = 'analog';
            Schan   = params.smartlab.tdt.co2(1);% analog chan were CO2 was recorded
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs     = ecog.analog.fs;
            Slabel   = ecog.analog.label{Schan};

            if Sdelay ~= 0 && correctdelay
                delaysamples = round(Sdelay*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
                % check size of vector
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end
                Ssignal = Ssignal2; clear Ssignal2
            end

            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.co2 = single(Ssignal); % collect
        end

        %% ---  O2 -----
        if isfield (params.smartlab.tdt,'o2') && ~isempty(params.smartlab.tdt.o2)
            Sdelay  = params.smartlab.tdt.o2(2);
            Stype   = 'analog';
            Schan   = params.smartlab.tdt.o2(1);% analog chan were O2 was recorded
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs     = ecog.analog.fs;
            Slabel  = ecog.analog.label{Schan};

            if  Sdelay ~= 0 && correctdelay
                delaysamples = round(Sdelay*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
                % check size of vector
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end
                Ssignal = Ssignal2; clear Ssignal2
            end

            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.o2 = single(Ssignal); % collect
        end


        %% ---  micro -----
        % micro can reveal important info
        if isfield (params.smartlab.tdt,'micro') && ~isempty(params.smartlab.tdt.o2)
            Sdelay  = params.smartlab.tdt.o2(2);
            Stype   = 'analog';
            Schan   = params.smartlab.tdt.o2(1);% analog chan were O2 was recorded
            Ssignal = ecog.analog.trial{1}(Schan,:);
            Sfs     = ecog.analog.fs;
            Slabel  = ecog.analog.label{Schan};

            if  Sdelay ~= 0 && correctdelay
                delaysamples = round(Sdelay*Sfs);
                signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
                signalpatch = fliplr(signalpatch);
                Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
                % check size of vector
                if numel(Ssignal2)>Ssignal
                    Ssignal2 = Ssignal2(1:numel(Ssignal));
                end
                Ssignal = Ssignal2; clear Ssignal2
            end

            if downsampling && fs~=Sfs
                newfs = fs; oldfs = Sfs;
                [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
                clear newfs oldfs
            end

            ecog.smartlab.tdt.o2 = single(Ssignal); % collect
        end

        % --- add extra info ---
        ecog.smartlab.tdt.time = (1:length(Ssignal))*1/fs ; % add time vector
        if downsampling && fs~=Sfs; ecog.smartlab.tdt.fs = fs; end % correct FS
    

    %% -------------------- end SMARTLAB to TDT linking  --------------------------------------------------

    %% get HR peak times ------
    % if exist('emg','var') && ~isempty(emg)
    %get the peaks
    %    [hrtimestamps] = BU_cleanEKG(emg,fs);%% use hrF if signal very noisy;use hr if signal ok
    %end

    %% >>>>> filtering & downsampling HR >>>>>
    if isfield(ecog,'osat')
        if isfield(ecog.osat, 'trial')
            osat = ecog.osat.trial;%%ecog.hr.trial{:};
            osat_fs = ecog.osat.fsample;
        end %plot(resp(1,1:20000))
        downsamplingS = 0;
        if downsamplingS && osat_fs ~= fs
            dataf = 500;%new fs
            fs_data = osat_fs;%old fs
            [upsamp_data,downsamp_data] = rat(dataf/fs_data,1e-10);
            osat = resample(double(osat)', upsamp_data, downsamp_data)'; % resample ECoG to dataf using
            clear dataf upsamp_data downsamp_data fs_data wn n a b
        end
        if length(osat) > length(eeg)
            disp("check lenghts");pause
            osat = osat(:,1:legth(eeg));
        end
    end

    %% >>>>> filtering & downsampling spo >>>>>
    if isfield(ecog,'spo')
        if isfield(ecog.spo, 'blood')
            spo      = ecog.spo.blood;%%ecog.hr.trial{:};
            spo_hr1  = ecog.spo.HRnatus;
            spo_hr2  = ecog.spo.HRsine;
            spo_fs   = ecog.spo.fsample;
        end %plot(resp(1,1:20000))
        downsamplingS = 0;
        if downsamplingS && spo_fs ~= fs
            dataf = 500;%new fs
            fs_data = spo_fs;%old fs
            [upsamp_data,downsamp_data] = rat(dataf/fs_data,1e-10);
            spo = resample(double(spo)', upsamp_data, downsamp_data)'; % resample ECoG to dataf using
            spo_hr1 = resample(double(spo_hr1)', upsamp_data, downsamp_data)'; % resample ECoG to dataf using
            spo_hr2 = resample(double(spo_hr2)', upsamp_data, downsamp_data)'; % resample ECoG to dataf using
            spo_fs  = dataf;
            clear dataf upsamp_data downsamp_data fs_data wn n a b
        end
        if length(spo) > length(eeg)
            disp("check lenghts");pause
            spo = spo(:,1:legth(eeg));
            spo_hr1 = spo_hr1(:,1:legth(eeg));
            spo_hr2 = spo_hr2(:,1:legth(eeg));
        end
    end

    %% ------- collect -------
    ecog1=[];
    %ecog.sampleON        = sampleON;
    %ecog.sampleOFF       = sampleOFF;
    %ecog1.eeg             = single(eeg);
    %ecog1.encodes         = encodes;
    %ecog1.encodesplus     = encodesplus;
    %ecog1.timestamps      = timestamps;
    %ecog1.labels_FSurf    = labels_FSurf;
    %ecog1.labels          = labels;
    %ecog1.fs              = fs;
    ecog1.params          = params;
    if downsampling && isfield(params,'sampleON') && ~isempty(params.sampleON)
        ecog1.params.sampleON  = params.sampleON;
        ecog1.params.sampleOFF = params.sampleOFF;
    else
        ecog1.params.sampleON  = 1;
        ecog1.params.sampleOFF =length(eeg);
    end

    ecog1.ftrip.trial{1}     = single(eeg);
    ecog1.ftrip.label_FSurf  = labels_FSurf;
    ecog1.ftrip.label        = labels;
    ecog1.ftrip.fsample      = fs;
    ecog1.ftrip.time{1}      = [0:length(eeg)-1]*[1/fs];
    ecog1.ftrip.nChans       = size(eeg,1);

    if isfield(ecog,'hr')
        ecog1.hr.trial      = single(ecog.hr.trial);
        ecog1.hr.fsample     = fs;
        ecog1.hr.label       = 'hr';
    end

    if isfield(ecog,'resp') && ~isempty(ecog.resp)
        ecog1.resp.trial{1} = resp;
        ecog1.resp.label    = params.resp.label; %params.resp_labels;
        ecog1.resp.fsample  = fs;
    end

    if isfield(ecog,'emg') && ~isempty(ecog.emg)
        ecog1.emg.trial{1}    = emg;
        if isfield(params.emg, 'chanlab'); ecog1.emg.label = params.emg.chanlab; end
        if isfield(params.emg, 'label'); ecog1.emg.label = params.emg.label; end
        ecog1.emg.fsample     = fs;
    end

    if isfield(ecog,'scalp')
        ecog1.scalp.trial{1} = scalp;
        ecog1.scalp.label    = ecog.scalp.label;
        ecog1.scalp.fsample  = fs;
        %ecog.timestamps.hr = timestamps_HR;
    end

    if isfield(ecog,'cranial') && ~isempty(ecog.cranial)
        ecog1.cranial.trial{1} = cranial;
        ecog1.cranial.label    = cranial_label';
        ecog1.cranial.fsample  = fs;
    end

    if isfield(ecog,'osat')
        ecog1.osat.trial{1} = osat;
        ecog1.osat.fsample  = fs;
    end

    if isfield(ecog,'spo')
        ecog1.spo.trial{1,1}     = spo;
        ecog1.spo.trial{1,2}     = spo_hr1;
        ecog1.spo.trial{1,3}     = spo_hr2;
        ecog1.spo.fsample        = spo_fs;
        ecog1.spo.encodes        = params.spo.encodes ;
    end

    if exist('timestamps') && ~isempty(timestamps)
        ecog1.timestamps.stim = timestamps;
        ecog1.timestamps.stimencodes = timestamps_encodes; %encodes;
        ecog1.timestamps.behaviour = timestamps_encodes; %encodes;
    end

    % sam added for gsm mixer output carbon dioxide flow percent
    % might have to add chanlabels similar to emg and resp if we are
    % getting more outputs from the gsm mixer 
    if isfield(ecog,'gsm_mixer') && ~isempty(ecog.gsm_mixer)
        ecog1.gsm_mixer.trial{1}    = gsm_mixer_ana;
        ecog1.gsm_mixer.label       = gsm_mixer_label;
        ecog1.gsm_mixer.fsample     = fs;
    end

    if isfield(params, 'SZ') && params.SZ==1
        %% visualize szr b4 saving - adjust params file
        ecognow = ecog1;
        cfg =[];
        cfg.channel = params.SZspike_chans';%{'LDa*'};%'all';%%only plot szr chans
        cfg.channel = params.SZonset_chans';%
        %%prepare extra chans to be added
        addekg =1;
        addresp =1;
        if addekg && addresp
            addthis=[hrF;resp];
            addthis_label=[ecognow.hr.label; ecognow.resp.label];
        end

        temp = {single(addthis)};
        ecognow.ftrip.trial{1} = [ecog1.ftrip.trial{1} ; addthis];
        ecognow.ftrip.nChans = ecog1.ftrip.nChans+size(addthis,1);
        ecognow.ftrip.label  =[ecog1.ftrip.label; addthis_label];
        cfg.channel = [cfg.channel ;addthis_label]
        %%discard bad chans
        ecognow.bad_chans = {};%{'LPc2'};
        params.SZtimestamps
        params.SZtimestamps_encodes
        %art=ecog_databrowser(ecognow,cfg)

        %% import neurologist notes (NEEDS WORK!)FOR NOW WRITE MANUALLY TO PARAMS
        %         filename = [params.filename '.txt'];%%'myfile01.txt';
        %         [A,delimiterOut] = importdata(filename);
        %         open A
        %         disp('NEUROLOGIST NOTES: make sure these strings are present: EEG ONSET CLINICAL ONSET,EEG OFFSET,CLINICAL OFFSET,EEG/CLINCAL OFFSET')
        %         dbstack;pause
        %         strings2find = {'EEG ONSET';'CLINICAL ONSET';'EEG/CLINCAL ONSET';'EEG OFFSET';'CLINICAL OFFSET';'EEG/CLINCAL OFFSET'};
        %          for i = 1:length(A)
        %             temp = A{i}
        %             for ii = 1:length(strings2find)
        %                X = strfind(temp,[strings2find{ii}]);
        %                if ~isempty
        %
        %             end
        %         end

        %%IMPORT directly from params for now
        %% save szr info
        ecog1.szr.timestamps         = params.SZtimestamps; % start & end of SZ
        ecog1.szr.onset_chans        = params.SZonset_chans; % contacts with clear sz onset activity (Ictal)
        ecog1.szr.spike_chans        = params.SZspike_chans; % contacts involved in sz (Ictal/Interictal) activity
        ecog1.szr.electropgraphic    = params.SZelectrographic;
        ecog1.szr.clinical           = params.SZclinicalsigns;
        ecog1.szr.timestamps_encodes = params.SZtimestamps_encodes;
    end

    %% start & stop of rec (utctime)
    if isfield(ecog, 'paramsTDT')
        ecog1.params.TDTrecON = ecog.paramsTDT;
    end

    if isfield(ecog, 'smartlab')
        ecog1.smartlab = ecog.smartlab;
    end

    ecog = ecog1; clear ecog1

    %% ------- add analog chans here --------
    clearvars -except ecog params

    save([params.directoryOUT params.filename '_ecog_epoch.mat'],'ecog');
    end  
    
end

%% ----------------------
%% pre-process SL data (uncomment and run b4 step 4)

% % % STEP 3.5
% [1] make SLtable
file_path = [params.smartlab.dir(1:end-3) 'mat'];%'/Volumes/F_BBio/RRST/patients/NS192/rawdata/logs/ori/SmartLab_5-24-2023_B14RRST.csv';
% % JA: now, go through JA_SLcsv2mat (adapted version of BU_SLcsv2mat, but actually
% % does mat2mat), and change column names
SLtable = JA_SLcsv2mat(file_path, ['/' params.filename '_SLtable']);% [SLtable] = BU_SLcsv2mat(file_path);
BU_InterpolateMissingSLtable(); % (JA: OPEN THIS) manually save SLtable to folder proc(eg. /Volumes/F_BBio/RRST/patients/NS195/rawdata/logs/process/B11_RRST_SLtable.m') 
% % then:Step4
% % [2] get timestamps for sink pulses (SLtable)
% figure;
SLtabledir = '/home/david/Documents/RRET/patients/NS221/vccp/logs/sl/NS221_B14_mouthVCCL_SLtable.mat';
columnX = 'DigInput1';%'Input2(tdt)'; %'Input_breathswitch' or 'Input1(ring)'... specify column name with input TTLs to get timestamps
[timeonsets,timeoffsets] = BU_getTimestampsfromSLtable(SLtabledir,columnX); % check timestamps
% %BU_getTimestampsfromSLtable_gpt % JA: idk what this is
% %manually save SL timestamps to excel 

%% ----------------------------------------------------------------------------
if Step4

    clear Step*

    %% load neural data (tdt block)

    load([params.directoryOUT params.filename '_ecog_epoch.mat']);

    %% replace params in case not updated

    ecog.params = params;

    %% COMPUTE TIMESTAMPS TO SINK AMPS
    % use digital triggers (sink-pulses sent to both TDT & SL if available)
    % all other data will be aligned with respect to this!!!!

    % compute sink times if not done yet
    disp('COMPUTE SINK TIMES IF NOT DONE ALREADY (ex. Experiment_log_NS196 -> B48_sink')

    if isfield(params.smartlab,'sinktimes') && isempty(params.smartlab.sinktimes)

        % get timestamps from SINK chan - TDT
        sinkpulsesTDT = ecog.smartlab.tdt.sinkpulses';
        
        % compute timestamps from SINK chan - SL
        [sinkpulsesSL, Input1] = BU_getTimestampsfromSLtable (params.smartlab.dir, params.smartlab.sinkID);

        if length(sinkpulsesSL)== length(sinkpulsesTDT)
            disp('GREAT!')
            disp('COPY TIMES MANUALLY TO SL_time and TDT_time IN EXCELL')
        else
            disp('NUMBER OF SINK PULSES IN SL & TDT DO NOT MATCH ---> CUTOFF SLtable OR TDT data so that both start at ~ same time')
        end

    end

    % load sink times if already done
    disp('LOAD SINK TIMES FROM EXCELL FILE (ex. Experiment_log_NS196 -> B48_sink')

    if isfield(params.smartlab,'sinktimes') && ~isempty(params.smartlab.sinktimes)

        % load sink timestamps from excel file (in secs)
        idx = 1:strfind(params.smartlab.sinktimes, '.xlsx');
        filenamenow = [params.smartlab.sinktimes(idx) 'xlsx']; % excell fullname including sheet
        sheet = params.smartlab.sinktimes(idx(end)+5:end);
        SinkTable = readtable(filenamenow, 'Sheet', 'NS221_B14_mouthVCCL_sink'); % load excell sheet with sink timestamps (in secs): 'SL_time', 'TDT_time'

        % provide titles if not picked up
        % SinkTable.Properties.VariableNames(1) = "SL_time";
        % SinkTable.Properties.VariableNames(2) = "TDT_time";
        % SinkTable.Properties.VariableNames(3) = "ET_time";
        % SinkTable.Properties.VariableNames(4) = "IDs";
        % SinkTable.Properties.VariableNames(5) = "Trial";
        % SinkTable.Properties.VariableNames(6) = "Pulse_num";

        % extract timestamps
        Cnames = SinkTable.Properties.VariableNames'; % 'SL_time', 'TDT_time'
        sinkpulsesTDT = SinkTable.TDT_time;
        sinkpulsesSL  = SinkTable.SL_time;

    end

    % % if isempty (sinkpulsesSL) && exist('Input1','var') && ~isempty(Input1)
    % %
    % %     if isfield(params.smartlab,'sinktimes') && isempty(params.smartlab.sinktimes)
    % %
    % %         sinkpulsesSL = [];
    % %         figure;subplot(311)
    % %         plot(Input1,'k');hold on;ylim([0 1.5]);
    % %         temp = [0 ;diff(Input1)];
    % %         [PKS,LOCS] = findpeaks(temp,"Threshold",.1,"MinPeakDistance",100);
    % %         vline(LOCS,'r:');
    % %         idx = find(diff(LOCS)>100*1.5);
    % %         idx = [1;idx+1];
    % %         % add numbers next to each vertical ttl bar
    % %         for i = 1:length(idx)
    % %             x = LOCS(idx(i));
    % %             line([x x], ylim, 'Color', 'g','LineStyle','--');  % This creates a vertical line at time x
    % %             text(x, 0.99, num2str(i), 'HorizontalAlignment', 'center');  % This adds the trial number next to the vertical line
    % %         end
    % %         %vline(LOCS(idx),'g');
    % %
    % %         % collect pulse timings (SL)
    % %
    % %         disp('Copy SL timestamps to excell sheet - collumn SL_time')
    % %         sinkpulsesSL = Time(LOCS(idx));
    % %
    % %         % get ttls sink pulses from TDT
    % %
    % %         sinkpulsesTDT = ecog.smartlab.tdt.sinkpulses';
    % %
    % %         % check pulses from SL & TDT do match
    % %
    % %         if numel(sinkpulsesSL) ~= numel(sinkpulsesTDT)
    % %             disp('SINK PULSES DO NOT MATCH -CHECK!!!');pause
    % %             subplot(312)
    % %             plot(Time,Input1,'k');hold on;ylim([0 1.5]);
    % %             for i = 1:length(sinkpulsesSL)
    % %                 x = sinkpulsesSL(i);
    % %                 line([x x], ylim, 'Color', 'r','LineStyle','--');  % This creates a vertical line at time x
    % %                 text(x, 0.99, num2str(i), 'HorizontalAlignment', 'center');  % This adds the trial number next to the vertical line
    % %             end
    % %         end
    % %
    % %         % adjust ttls from both amps to SAME CLOCK
    % %         subplot(313)
    % %         % adjust SL clock
    % %         diffnow = sinkpulsesTDT(1) - sinkpulsesSL(1);
    % %         tempSLpulses = sinkpulsesSL + diffnow;
    % %         tempSLpulses(1)
    % %         vline(tempSLpulses,'k');
    % %         hold on;ylim([0 1.5]);
    % %
    % %         % overlie TDT clock
    % %         vline(sinkpulsesTDT,'r:');
    % %
    % %     else
    % %
    % %         % load sink timestamps from excel file (in secs)
    % %         idx = 1:strfind(params.smartlab.sinktimes, '.xlsx');
    % %         filenamenow = [params.smartlab.sinktimes(idx) 'xlsx']; % excell fullname including sheet
    % %         sheet = params.smartlab.sinktimes(idx(end)+5:end);
    % %         SinkTable = readtable(filenamenow, 'Sheet', sheet); % load excell sheet with sink timestamps (in secs): 'SL_time', 'TDT_time'
    % %
    % %         % extract timestamps
    % %         Cnames = SinkTable.Properties.VariableNames'; % 'SL_time', 'TDT_time'
    % %         sinkpulsesTDT = SinkTable.TDT_time;
    % %         sinkpulsesSL  = SinkTable.SL_time;
    % %     end
    % % end

    %% EPOCH TDT DATA BASED ON TIMESTAMPS ABOVE
    % 'ringexp1' ttl will be time 0s in every epoch

    % epoching NEURAL data first
    ecog.epochs =[];

    % sam added for vccp protocol
    if strcmp(params.task,'VCCP')
        %pre = 10; pos = 130; % select data around event ('trialstart' which is time=0s ) % consdiering 140 as 120 seconds after the last trial start but still
        pre = 10; pos = 130;
        % there was an error in RRBO if i had no pre seconds, so taking 10
        % seconds since every breath has a 8 second context required
        numringexp    = 0; % no ring
        sinkID = 'trialstart'; %only TTL for start and stop of trial
        stopID = 'trialstop';

    elseif strcmp(params.task,'RRET')
        pre = 60; pos = 160; % *** change back to 160s % select time window around event - 'ringexp1' is time 0s
        numringexp = 20; % ring expands 10 times [***JA: different # of expansions on different trials... doesn't seem to matter]
        sinkID = 'ringexp1';
    end


    index = find(strcmp(SinkTable.IDs, sinkID)==1); % select event type (IDs)
    stop_index = find(strcmp(SinkTable.IDs, stopID)==1);
    timestampstdt = sinkpulsesTDT(index);
    timestampsSL  = sinkpulsesSL(index);
    timestampstdt_stop = sinkpulsesTDT(stop_index);
    timestampsSL_stop = sinkpulsesSL(stop_index);

    TDT_trial_dur = timestampstdt_stop - timestampstdt;

    % epoch eeg
    vector1 = ecog.ftrip.trial{1};
    fs1 = ecog.ftrip.fsample;
    [datanow,timenow,timetracker] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
    %%%[datanow,timenow] = BU_epochSignals(vector1, fs1, timestampstdt, pre, pos);

    % collect epoched eeg
    ecog.epochs.eeg = [];
    ecog.epochs.eeg  = datanow;
    ecog.epochs.time = timenow;
    ecog.epochs.timetracker = timetracker;
    ecog.epochs.fs = fs1;
    ecog.epochs.eeglabels = ecog.ftrip.label;
    ecog.epochs.eeglabels_FSurf = ecog.ftrip.label_FSurf;
    ecog.epochs.sinkID = sinkID; %(ex.'ringexp1' or 'breathswithc));

 
    %% check types of events  
    IDsall = unique (SinkTable.IDs);

    if strcmp(params.task,'RRST') || strcmp(params.task,'RRET')

        containsKeypress = ismember('keypress', IDsall);
        if containsKeypress
            indexref = find(strcmp(SinkTable.IDs, 'keypress')==1); % keypress
            timenow = mean(sinkpulsesTDT(indexref)-tibeginsmesptampstdt); % calculate the mean of all trials
            idx = findnearest (ecog.epochs.time,timenow); % find closest point in time
            ecog.epochs.time_keypress = ecog.epochs.time(idx);
            ecog.epochs.time_keypressTrials = sinkpulsesTDT(indexref)-timestampstdt;
        end

        containsGetready = ismember('getready', IDsall);
        if containsGetready
            indexref = find(strcmp(SinkTable.IDs, 'getready')==1); % get ready to inhale (ring changes color)
            timenow = mean(sinkpulsesTDT(indexref)-timestampstdt); % calculate the mean of all trials
            idx = findnearest (ecog.epochs.time,timenow); % find closest point in time
            ecog.epochs.time_getready = ecog.epochs.time(idx);
            ecog.epochs.time_getreadyTrials = sinkpulsesTDT(indexref)-timestampstdt;
        end

        % JA commented below: unnecessary?
        % for i = 2:numringexp
        %     ringexp_name = ['ringexp', num2str(i)];
        %     indexref = find(strcmp(SinkTable.IDs, ringexp_name)==1);
        %     %indexref = [indexref; index];  % Append the found index to the array
        %     timenow = mean(sinkpulsesTDT(indexref)-timestampstdt); % calculate the mean of all trials
        %     idx = findnearest (ecog.epochs.time,timenow); % find closest point in time
        %     ecog.epochs.(['time_', ringexp_name]) = ecog.epochs.time(idx); % collect
        %     ecog.epochs.(['time_', ringexp_name,'Trials']) = sinkpulsesTDT(indexref)-timestampstdt; % collect
        % end
    end

    if isfield(ecog,'resp') && ~isempty(ecog.resp)
        ecog.epochs.resp =[];
        vector1 = ecog.resp.trial{1}(1,:); %samruddhi changed from trial{1} tp trial{1}(1,:) to choose one resp channel/can be done in earlier steps but had that commented
        fs1 = ecog.resp.fsample;
      
        vector1 = vector1*-1;   % make sure the inh is downwards 
        maxAbsVal = max(abs(vector1));% scaled from +1:1
        vector1 = vector1 / maxAbsVal;

        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.resp  = datanow;
        ecog.epochs.resplabels = ecog.resp.label;
    end

    %% Sam added for additional analogs from TDT

    if isfield(ecog,'gsm_mixer') && ~isempty(ecog.gsm_mixer)
        ecog.epochs.gsm_mixer =[];
        vector1 = ecog.gsm_mixer.trial{1}; 
        fs1 = ecog.gsm_mixer.fsample;
        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.gsm_mixer  = datanow;
        ecog.epochs.gsmmixlabel = ecog.gsm_mixer.label;
    end
   %%
    if isfield(ecog,'cranial') && ~isempty(ecog.cranial)
        ecog.epochs.cranial =[];
        vector1 = ecog.cranial.trial{1};
        fs1 = ecog.cranial.fsample;
        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.cranial = datanow;
        ecog.epochs.craniallabels = ecog.cranial.label;
    end

    if isfield(ecog,'emg') && ~isempty(ecog.emg)
        ecog.epochs.emg =[];
        vector1 = ecog.emg.trial{1};
        fs1 = ecog.emg.fsample;
        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.emg = datanow;
        ecog.epochs.emglabels = ecog.emg.label;
    end

    % EPOCH INHALATION DURATIONS (INHONOFF) OR INH ONSETS OFFSETS
    % Pulses sent from SL-digout1 to TDT-anag chan
    % If the correct TDT program was used it will digitize output1 analog signal - sending times
    % in secs) for inh onsets/offsets but this NOT working
    % great

    if isfield(ecog.smartlab.tdt,'inhON')
        if isempty(ecog.smartlab.tdt.inhON)
            %% NEEDS WORK
        end
    end

    if isfield(ecog.smartlab.tdt,'inhON') && ...
            ~isempty(ecog.smartlab.tdt.inhON)

        % epoch inh onsets
        ecog.epochs.inhonsets_tdt = cell(length(timestampstdt), 1);
        timestamps = ecog.smartlab.tdt.inhON;
        for i = 1:length(timestampstdt)
            idx = find(timestamps >= timestampstdt(i)-pre & timestamps <= timestampstdt(i)+pos);
            temp = timestamps(idx) - timestampstdt(i); % time 0 is ring1exp
            % correct delay
            temp = temp-params.smartlab.tdt.inhON(2);%SLoutput also has .2+.04 delay
            ecog.epochs.inhonsets_tdt{i} =temp;
        end

        % epoch inh offsets
        if ~isempty(ecog.smartlab.tdt.inhOFF)
            ecog.epochs.inhoffsets_tdt = cell(length(timestampstdt), 1);
            timestamps = ecog.smartlab.tdt.inhOFF;
            for i = 1:length(timestampstdt)
                idx = find(timestamps >= timestampstdt(i)-pre & timestamps <= timestampstdt(i)+pos);
                temp = timestamps(idx) - timestampstdt(i); % time 0 is ring1exp

                % correct delay
                temp = temp - params.smartlab.tdt.inhONOFF(2);%SLoutput has a .2+.04 delay
                ecog.epochs.inhoffsets_tdt{i} = temp;
            end
        end
    end


    %% SL-TDT analogs (pressure)

    if isfield(ecog.smartlab.tdt,'pressure') && ~isempty(ecog.smartlab.tdt.pressure)
        ecog.epochs.pressure_tdt =[];
        vector1 = ecog.smartlab.tdt.pressure;
        fs1 = 500;%ecog.smartlab.tdt.fs;

        % corrrect delay (this was done on step 3 ... check again)

        %idx1 = round(params.smartlab.tdt.pressure(2)*fs1+1) : numel(vector1); % shift data forward idx
        %idx2 = numel(vector1)-idx1(1)+2:numel(vector1); % fill up end to make it same size window
        %idx2 = sort(idx2,"descend");
        %vector1 =  vector1([idx1 idx2]);

        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.pressure_tdt = datanow;
    end

    %% SL-TDT analogs (Volume)

    if isfield(ecog.smartlab.tdt,'volume') && ~isempty(ecog.smartlab.tdt.volume)
        ecog.epochs.volume_tdt =[];
        vector1 = ecog.smartlab.tdt.volume;
        fs1 = 500;%ecog.smartlab.tdt.fs;

        % corrrect delay (this was done on step 3... check again possible)

        idx1 = round(params.smartlab.tdt.volume(2)*fs1+1) : numel(vector1); % starting sample
        idx2 = numel(vector1)-idx1(1)+2:numel(vector1); % fill up end to make it same size window
        idx2 = sort(idx2,"descend");
        vector1 =  vector1([idx1 idx2]);

        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.volume_tdt = datanow;
    end

     %% SL-TDT analogs (Flow)

    if isfield(ecog.smartlab.tdt,'flow') && ~isempty(ecog.smartlab.tdt.flow)
        ecog.epochs.flow_tdt =[];
        vector1 = ecog.smartlab.tdt.flow;
        fs1 = 500;%ecog.smartlab.tdt.fs;

        % corrrect delay (this was done on step 3... check again possible)

        idx1 = round(params.smartlab.tdt.flow(2)*fs1+1) : numel(vector1); % starting sample
        idx2 = numel(vector1)-idx1(1)+2:numel(vector1); % fill up end to make it same size window
        idx2 = sort(idx2,"descend");
        vector1 =  vector1([idx1 idx2]);

        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.flow_tdt = datanow;
    end

    %% SL-TDT analogs (CO2)

    if isfield(ecog.smartlab.tdt,'co2') && ~isempty(ecog.smartlab.tdt.co2)
        ecog.epochs.co2_tdt =[];
        vector1 = ecog.smartlab.tdt.co2;
        fs1 = 500;%ecog.smartlab.tdt.fs;

        % corrrect delay (this was done on step 3... check again possible)
        % sam - this step has been already run in step 3 then why repeat?

        % delaynow = params.smartlab.tdt.co2(2); % windows delay + Instacall delay
        % idx1 = round(delaynow*fs1)+1 : numel(vector1); % starting sample
        % idx2 = numel(vector1)-idx1(1)+2:numel(vector1); % fill up end to make it same size window
        % idx2 = sort(idx2,"descend");
        % vector1 =  vector1([idx1 idx2]);

        %epoch
        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.co2_tdt = datanow;
    end

     %% SL-TDT analogs (O2)

    if isfield(ecog.smartlab.tdt,'o2') && ~isempty(ecog.smartlab.tdt.o2)
        ecog.epochs.o2_tdt =[];
        vector1 = ecog.smartlab.tdt.o2;
        fs1 = 500;%ecog.smartlab.tdt.fs;

        % corrrect delay (this was done on step 3... check again possible)

        delaynow = params.smartlab.tdt.o2(2); % windows delay + Instacall delay
        idx1 = round(delaynow*fs1)+1 : numel(vector1); % starting sample
        idx2 = numel(vector1)-idx1(1)+2:numel(vector1); % fill up end to make it same size window
        idx2 = sort(idx2,"descend");
        vector1 =  vector1([idx1 idx2]);

        %epoch
        [datanow] = BU_EpochDataTrials(vector1, fs1, fs1, timestampstdt, timestampstdt_stop);
        ecog.epochs.o2_tdt = datanow;
    end


    %% SMARTLAB DIRECT  (XML or CSV)
    % 1) smooth SL data
    % 2) resample (100hz > 500hz)
    % 2) epoch

    % SL - Pressure
    SLfs = params.smartlab.fs;

    % extract variables to workspace (usually 8 columns)
    if ~exist ("Pressure",'var')
        
        load([fileparts(params.smartlab.dir) '/' params.filename '_SLtable']);% JA commented old: load(params.smartlab.dir);% load data (SLdata)
        
        column_names = SLtable.Properties.VariableNames';% columns titles (Time,Flow,Volumen,Pressure,...)
        for j = 1:numel(column_names)
            column_name = column_names{j};
            column_values = SLtable.(column_name);

            longname = strfind(column_name,'(');
            if ~isempty(longname)
                column_name = column_name(1:longname-1);
            end

            eval([column_name ' = column_values;']);
        end
    end

    %% CHECK SL TIME VECTOR  
    % if not equal to TDT time vector then change timestamps to match
        
    if exist('Time')
        if Time(1) >= 1/SLfs
            response = questdlg('CAREFULL - TDT time vecotr does not match SL time vector.... CORRECT?', 'Offset Adjustment', 'Yes', 'No', 'Yes');
            if strcmp(response,'Yes')
               if Time(1)<0;timestampsSL = timestampsSL-Time(1);end %if negative add it on
               if Time(1)>0;timestampsSL = timestampsSL+Time(1);end %if positive substract it 
            end
        end
    end


    if exist ("Pressure",'var')

        figure;
        plot(Pressure,'b');hold on
        % apply a low-pass filter to the Pressure signal to smooth it a bit

        cutoff_freq = 3;  % cut-off frequency of the low-pass filter, in Hz
        order = 2;  % order of the filter
        [b, a] = butter(order, cutoff_freq/(SLfs/2), 'low');
        Pressure = filtfilt(b, a, Pressure);  % Zero-phase filtering
        clear cutoff_freq order b a
        plot(Pressure,'r');hold on

        %% Adjust offset if needed
        adjustoffset = 0;
        if adjustoffset
            % identify the Lowest Values
            percentage = 10;% lowest 5% of values (adjust the percentage as needed).
            numValues = numel(Pressure);
            if median(Pressure<0)
                sortedSignal = sort(Pressure, 'descend');% for negative signals use 'descend'
            end
            lowestValues = sortedSignal(1:round(numValues * percentage / 100));
            medianLowest = median(lowestValues);% Calculate the Median of the Largest Values
            adjustedPressure = Pressure - medianLowest;%Shift the Signal
            plot(adjustedPressure,'m');hline(0,'m')
            Pressure = adjustedPressure;% **
            plot(Pressure,'m');hold on
            disp('Pressure signal (magenta) LOOKING GOOD?');pause;close all hidden
        end
        
        %% upsample freq to 500hz and epoch

        oldfs = SLfs; newfs = ecog.epochs.fs;
        [epochsnow] = BU_EpochDataTrials(Pressure, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);%HERE THIS SHOULD GO BY TIMES IN SL (INPUT TIMEVECTOR
        ecog.epochs.pressure = epochsnow;
        ecog.epochs.pressurelabels = 'Pressure (cmH2O)';

        %% compute inh onsets/offsets based on SL-pressure signal

        ecog.epochs.inhonsets  = {};
        ecog.epochs.inhoffsets = {};

        %% use SL-Flow if present - else use SL-Pressure
        if isduration(Time)
        Time = seconds(Time);
        end

        [onsets, offsets] = BU_autofindInhOnsOffs(IFlow, Pressure, Time, 'trial', CO2, DigInput1, O2);
        %[onsets, offsets] = BU_detect_InhonsetsFromPressure(IFlow, Time, SLfs);
        %[onsets, offsets] = BU_detect_InhonsetsFromPressure(Flow*-1, Time, SLfs);
        %[onsets, offsets] = BU_detect_InhonsetsFromPressure(Pressure, Time, SLfs);
        %%[onsets, offsets] = BU_detect_InhonsetsFromPressure(epochsnow, ecog.epochs.time, SLfs);

        % JA *** if Time(1)~=0, the onsets/offsets are messed up... not
        % sure why this happened on NS206_B12

        % BREAK POINT TO EXAMINE AND NOTE BREATH NUMBERS AND GOOD/BAD

        if Time(1)~=0
            onsets = onsets-Time(1);
            offsets = offsets-Time(1);
        end

        % JA: manually select inh onsets+offsets for each trial (write down
        % conditions too for later in BXX_trialsinfo in exp log)
        % firstBreaths = [13 37 61 87 111 137 179 202 225 246 273 296 317 343 364 388 410 433]';
        % lastBreaths = [33 57 80 107 131 157 199 221 244 266 293 316 340 363 386 409 432 455]';
        % trialBreathIdx = [firstBreaths lastBreaths];
        % 
        % firstChal = [0 42 66 92 116 143 184 207 229 251 0 301 326 348 371 396 419 441];
        % firstReli = [0 53 76 103 127 152 194 216 239 262 0 311 335 359 381 405 430 451];

        firstBreaths = [129 177 216 233 266 294 324 344 360 385 408 432 456]';
        lastBreaths = [172 205 228 262 289 322 341 357 381 405 429 453 495]';
        trialBreathIdx = [firstBreaths lastBreaths];

        %% sam commented for vccp

        % firstChal = [0 59 85 114 137 160 190 212 233];
        % firstReli = [0 71 95 125 146 170 199 0 0];
        % 
        % %use this to compute number of each breath condition, then paste into trialsinfo spreadsheet
        % numTrials = length(epochsnow);
        % buddy = zeros(numTrials, 3);
        % for trialI = 1:numTrials
        %     % number of controls
        %     if firstChal(trialI) ~= 0
        %         buddy(trialI, 1) = firstChal(trialI) - trialBreathIdx(trialI, 1);
        %     else
        %     buddy(trialI, 1) = diff(trialBreathIdx(trialI, :))+1;
        %     end
        %     % number of challenge
        %     if firstChal(trialI) ~= 0
        %         if firstReli(trialI) ~= 0
        %             buddy(trialI, 2) = firstReli(trialI) - firstChal(trialI);
        %         else
        %             buddy(trialI, 2) = trialBreathIdx(trialI, 2) - firstChal(trialI)+ 1;
        %         end
        %     end
        %     %number of relief
        %     if firstReli(trialI) ~= 0
        %         buddy(trialI, 3) = trialBreathIdx(trialI, 2) - firstReli(trialI) +1;
        %     end
        % end
        % 
        % BREAK POINT TO COPY BUDDY VAR TO TRIALSINFO TABLE ON EXP LOG
%% sam commented above for vccp
        % epoch onsets/offsets based on threshold crossings
        for i = 1:numel(timestampsSL)
            % idx = find(onsets >= timestampsSL(i)-pre & onsets <=
            % timestampsSL(i)+pos); % JA replaced with manual selection
            idx = trialBreathIdx(i, 1):trialBreathIdx(i, 2);
            ecog.epochs.inhonsets{i,1} = onsets(idx)-timestampsSL(i);
            ecog.epochs.inhoffsets{i,1} = offsets(idx)-timestampsSL(i);
        end
    end

    %% SL DIRECT - DigOut1 - Inhalation duration
    % Inhalation onset is when pressure goes <-.1
    % Inhalation offset is when pressure goes <-.1
    % there's a 30-50ms delay

    if exist ("Output1",'var') && isfield(params.smartlab.tdt,'inhONOFF') && ~isempty(isfield(params.smartlab.tdt.inhONOFF))
        vector1 = Output1;

        % correct delay
        delaynow = params.smartlab.tdt.inhONOFF(2) - .2; % windows delay only (plugged directly not through usb -Instacal )
        idx1 = round(delaynow * SLfs) + 1 : numel(vector1); % starting idx
        idx2 = numel(vector1) - idx1(1) + 2 : numel(vector1); % fill up end to make it same size window
        idx2 = sort(idx2,"descend");
        vector1 =  vector1([idx1 idx2]);

        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;
        [epochsnow] = BU_EpochDataTrials(vector1, oldfs,  newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.inhdur = epochsnow;
        ecog.epochs.inhdurlabels = 'InhONOFF';

        %% analog to digital conversion
        % NEEDS WORK TO DIGITIZE THIS BUT USING METHOD ABOVE that computes
        % onsets from pressure signal itself
        %         temp = Output1;
        %         temp = [0;diff(temp) == 1];
        %         plot(time,Output1); hold on; plot(time,temp,'k')
        %         findpeaks(temp,time,'Threshold',.1,'MinPeakProminence',.1,'MinPeakDistance',ecog.epochs.fs)
        %         hold on;
        %         plot(time,Pressure)
        %        findpeaks(temp,time)
        %
        %        hold on;plot(time,Output1,'r')
    end

    %% SL DIRECT - FLOW

    if exist ("Flow",'var') || exist ("Airflow",'var')
       
        if exist ("Airflow",'var');Flow = Airflow;end
        if exist ("IFlow",'var');Flow = IFlow;end
        
      
        % apply a low-pass filter to the Flow signal to zero and smooth the signal

        cutoff_freq = 3;  % cut-off frequency of the low-pass filter, in Hz
        order = 2;  % order of the filter
        [b, a] = butter(order, cutoff_freq/(SLfs/2), 'low');
        FlowF = filtfilt(b, a, Flow);  % Zero-phase filtering
        clear cutoff_freq order b a
        figure; plot(Flow,'b');hold on;plot(FlowF,'r');hold on
        Flow = FlowF;clear FLowF;

        % Adjust offset if needed
        response = questdlg('Adjust offset if needed?', 'Offset Adjustment', 'Yes', 'No', 'Yes');
        if strcmp(response,'Yes')
            % identify the Lowest Values
            percentage = 5;% lowest 10% of values. Adjust the percentage as needed(5 or 10 usually)
            numValues = numel(Flow);
            sortedSignal = sort(Flow, 'descend');% for negative signals use 'descend'
            lowestValues = sortedSignal(1:round(numValues * percentage / 100));
            medianLowest = median(lowestValues);% Calculate the Median of the Largest Values
            adjustedFlow = Flow - medianLowest;%Shift the Signal
            plot(adjustedFlow,'m');
            Flow = adjustedFlow;% **
            disp('Flow signal (magenta) LOOKING GOOD?');pause
        end

        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(vector1, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.flow = epochsnow;
        ecog.epochs.flowlabels = 'iFlow (LPM)';
        %ecog.epochs.flowlabels = 'Flow (LPM)';
    end

    % JA split IFlow from above
    if  exist ("IFlow",'var')
        Flow = IFlow;
        % apply a low-pass filter to the Flow signal to zero and smooth the signal

        cutoff_freq = 3;  % cut-off frequency of the low-pass filter, in Hz
        order = 2;  % order of the filter
        [b, a] = butter(order, cutoff_freq/(SLfs/2), 'low');
        FlowF = filtfilt(b, a, Flow);  % Zero-phase filtering
        clear cutoff_freq order b a
        figure; plot(Flow,'b');hold on;plot(FlowF,'r');hold on
        Flow = FlowF;clear FLowF;

        % Adjust offset if needed
        response = questdlg('Adjust offset if needed?', 'Offset Adjustment', 'Yes', 'No', 'Yes');
        if strcmp(response,'Yes')
            % identify the Lowest Values
            percentage = 5;% lowest 10% of values. Adjust the percentage as needed(5 or 10 usually)
            numValues = numel(Flow);
            sortedSignal = sort(Flow, 'descend');% for negative signals use 'descend'
            lowestValues = sortedSignal(1:round(numValues * percentage / 100));
            medianLowest = median(lowestValues);% Calculate the Median of the Largest Values
            adjustedFlow = Flow - medianLowest;%Shift the Signal
            plot(adjustedFlow,'m');
            Flow = adjustedFlow;% **
            disp('Flow signal (magenta) LOOKING GOOD?');pause
        end

        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(Flow, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.iflow = epochsnow;
        ecog.epochs.iflowlabels = 'IFlow (LPM)';
        %ecog.epochs.flowlabels = 'Flow (LPM)';
    end
    
    % JA added:
    if  exist ("EFlow",'var')
        Flow = EFlow;
        % apply a low-pass filter to the Flow signal to zero and smooth the signal

        cutoff_freq = 3;  % cut-off frequency of the low-pass filter, in Hz
        order = 2;  % order of the filter
        [b, a] = butter(order, cutoff_freq/(SLfs/2), 'low');
        FlowF = filtfilt(b, a, Flow);  % Zero-phase filtering
        clear cutoff_freq order b a
        figure; plot(Flow,'b');hold on;plot(FlowF,'r');hold on
        Flow = FlowF;clear FLowF;

        % Adjust offset if needed
        response = questdlg('Adjust offset if needed?', 'Offset Adjustment', 'Yes', 'No', 'Yes');
        if strcmp(response,'Yes')
            % identify the Lowest Values
            percentage = 5;% lowest 10% of values. Adjust the percentage as needed(5 or 10 usually)
            numValues = numel(Flow);
            sortedSignal = sort(Flow, 'descend');% for negative signals use 'descend'
            lowestValues = sortedSignal(1:round(numValues * percentage / 100));
            medianLowest = median(lowestValues);% Calculate the Median of the Largest Values
            adjustedFlow = Flow - medianLowest;%Shift the Signal
            plot(adjustedFlow,'m');
            Flow = adjustedFlow;% **
            disp('Flow signal (magenta) LOOKING GOOD?');pause
        end

        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(Flow, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.eflow = epochsnow;
        ecog.epochs.eflowlabels = 'EFlow (LPM)';
        %ecog.epochs.flowlabels = 'Flow (LPM)';
    end

    %% SL DIRECT - Volume

    if exist ("Volume",'var')
        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;
        [epochsnow] = BU_EpochDataTrials(vector1, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.volume = epochsnow;
        ecog.epochs.volumelabels = 'Volume (mL)';
    end

    %% SL DIRECT - CO2

    if exist ("CO2",'var')

        vector1 = CO2;

        % correct delay
        delaynow = params.smartlab.co2delay; % within SL delay (due to CO2 sensor delay)
        idx1 = round(delaynow * SLfs) + 1 : numel(vector1); % starting idx
        idx2 = numel(vector1) - idx1(1) + 2 : numel(vector1); % fill up end to make it same size window
        idx2 = sort(idx2,"descend");
        vector1 =  vector1([idx1 idx2]);

        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;
        [epochsnow] = BU_EpochDataTrials(vector1, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.co2 = epochsnow;
        ecog.epochs.co2labels = 'CO2 concentration (%)';

        %[epochs_amp1, epochs_amp2_resampled] = BU_SyncAmps(amp1_data, amp2_data, timestamps1, timestamps2, fs1, fs2)
        %[vector1_epochs,vector2_epochs,time_vector,vector2plus_epochs,TDTtimestamps] = ...
        %BU_AlignDatafrom2amps_timestamps (vector1,vector2,vector2plus,sampling_rate1,sampling_rate2,alignMethod)
    end

     %% SL DIRECT - O2
    if exist ("O2",'var')

        vector1 = O2;

        % correct delay
        delaynow = params.smartlab.o2delay; % windows/sensor delay only (no Instacall delay)
        idx1 = round(delaynow * SLfs) + 1 : numel(vector1); % starting idx
        idx2 = numel(vector1) - idx1(1) + 2 : numel(vector1); % fill up end to make it same size window
        idx2 = sort(idx2,"descend");
        vector1 =  vector1([idx1 idx2]);

        % correct delay
        oldfs = SLfs; newfs = ecog.epochs.fs;
        [epochsnow] = BU_EpochDataTrials(vector1, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.o2 = epochsnow;
        ecog.epochs.o2labels = 'O2 concentration (%)';
    end

    %% JA added: DigInput1

    if  exist ("DigInput1",'var')
        diginput = DigInput1;
        
        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(diginput, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.diginput1 = epochsnow;
        ecog.epochs.diginput1labels = 'DigInput1';
    end

    %% JA added: DigInput2

    if  exist ("DigInput2",'var')
        diginput = DigInput2;
        
        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(diginput, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.diginput2 = epochsnow;
        ecog.epochs.diginput2labels = 'DigInput2';
    end

    %% JA added: SpO2
    
    if  exist ("SpO2",'var')
        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(SpO2, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.spo2 = epochsnow;
        ecog.epochs.spo2labels = 'SpO2';
    end

    %% JA added: HeartRate
    
    if  exist ("HeartRate",'var')
        
        % upsample & epoch
        oldfs = SLfs; newfs = ecog.epochs.fs;

        [epochsnow] = BU_EpochDataTrials(HeartRate, oldfs, newfs, timestampsSL, timestampsSL_stop, TDT_trial_dur);
        ecog.epochs.heartrate = epochsnow;
        ecog.epochs.heartratelabels = 'HeartRate';
    end
    
    %% JA added: plot all epoched vars in all trials to check epoching
    
    ntrs = length(ecog.epochs.eeg);
    T = ecog.epochs.time;
    allFields = fieldnames(ecog.epochs); 
    epochFields = allFields(~ismember(allFields, {'inhonsets','inhoffsets'}));
    %epochFields = fieldnames(ecog.epochs);
    for trialI = 1:ntrs
        nexttile
        isEpochedVar = zeros(size(epochFields));
        for fieldI = 1:length(epochFields)
            buddy = ecog.epochs.(epochFields{fieldI});
            if iscell(buddy) && (min(size(buddy)==[ntrs 1])==1) && length(buddy{trialI}) == length(T) && ~strcmp(epochFields{fieldI}, 'eeg')
                plot(T{trialI,1}, buddy{trialI})
                hold on
                isEpochedVar(fieldI) = 1;
            end
        end
        xline(ecog.epochs.inhonsets{trialI})
        legend(epochFields(find(isEpochedVar)))
    end

    

    %%
    % %     if exist('breathswitch','var') && ~isempty(breathswitch)
    % %         vector1 = double(ecog.smartlab.tdt.breathswitch); % neural amp (500 Hz sampled data)
    % %         vector2 = double(breathswitch'); % Smartlab (10 Hz sampled data)
    % %         fs1 = ecog.smartlab.tdt.fs;
    % %         fs2 = 100; % params.smartlab.fs;
    % %         vectorplus(1,:) = double(balonONOFF'); % smartlab (10 Hz sampled data)
    % %         vectorplus(2,:) = double(Pressure'); % smartlab (10 Hz sampled data)
    % %         vectorplus(3,:) = double(CO2'); % smartlab (10 Hz sampled data)
    % %         vectorplus(4,:) = double(O2'); % smartlab (10 Hz sampled data)
    % %
    % %         [vector1_epochs,vector2_epochs,vector_time,vectorplus_epochs] = ...
    % %             BU_AlignDatafrom2amps_timestamps(vector1,vector2,vectorplus,fs1,fs2,'digital');
    % %     end

    %% if no digital input triggers detected align with digital output (balonONOFF)

    % %     if exist('balonONOFF','var') && ~isempty(balonONOFF) && ~exist('vector1_epochs','var')
    % %         vector1 = double(ecog.smartlab.tdt.balonONOFF); % neural amp (500 Hz sampled data)
    % %         vector2 = double(balonONOFF'); % smartlab (10 Hz sampled data)
    % %         fs1 = ecog.smartlab.tdt.fs;
    % %         fs2 = 100; % params.smartlab.fs;
    % %         vectorplus(1,:) = Pressure; % smartlab (10 Hz sampled data)
    % %         vectorplus(2,:) = CO2; % smartlab (10 Hz sampled data)
    % %         vectorplus(3,:) = O2; % smartlab (10 Hz sampled data)
    % %
    % %         [TDT_epochs,SL_epochs,vector_time,SL_epochsplus,TDTtimestamps] = ...
    % %             BU_AlignDatafrom2amps_timestamps (vector1,vector2,vectorplus,fs1,fs2,'digital');
    % %     end

    %% epoch tdt analog chans - final comparison
    % %     TDT_epochsplus =[];
    % %     TDTtimestamps;
    % %
    % %     %% ana1
    % %     ana1 = ecog.smartlab.tdt.breathswitch;
    % %     [ana1_epochs] = BU_epochSignal(ana1, fs1, TDTtimestamps, pre, pos);
    % %
    % %     plotting_tr = 1;
    % %     if plotting_tr
    % %         for i = 1:76;clf;
    % %             plot(vector_time,TDT_epochs(i,:));hold on;
    % %             plot(vector_time,ana1_epochs(i,:),'r');
    % %             legend('tdt balonONOFF', 'tdt-breathswitch');
    % %             title(['epoch=' num2str(i)]);
    % %             pause(2)
    % %         end
    % %     end
    % %
    % %     %%ana2
    % %     ana2 = ecog.smartlab.tdt.balon_ONOFF;
    % %     [ana2_epochs] = BU_epochSignal(ana2, fs1, TDTtimestamps, pre, pos);
    % %
    % %     if plotting_tr
    % %         for i = 1:76;clf;
    % %             plot(vector_time,TDT_epochs(i,:));hold on;
    % %             plot(vector_time,ana2_epochs(i,:),'r:');
    % %             legend('tdt balonONOFF', 'tdt balonONOFF');
    % %             title(['epoch #' num2str(i)]);
    % %             pause(2)
    % %         end
    % %     end
    % %
    % %     %%ana3
    % %     ana3 = ecog.smartlab.tdt.pressure;
    % %     [ana3_epochs] = BU_epochSignal(ana3, fs1, TDTtimestamps, pre, pos);
    % %     pressure = squeeze(SL_epochsplus(1,:,:));
    % %
    % %     if plotting_tr
    % %         for i = 1:76;clf;
    % %             plot(vector_time,ana3_epochs(i,:),'k');hold on;
    % %             plot(vector_time,pressure(i,:),'r:');
    % %             legend('TDT-pressure', 'SL-pressure');
    % %             title(['epoch #' num2str(i)]);
    % %             pause(2)
    % %         end
    % %     end
    % %
    % %
    % %     %%ana4
    % %     ana4 = ecog.smartlab.tdt.co2;
    % %     [ana4_epochs] = BU_epochSignal(ana4, fs1, TDTtimestamps, pre, pos);
    % %     co2 = squeeze(SL_epochsplus(2,:,:));

    % %     if plotting_tr
    % %         for i = 1:76;clf;
    % %             plot(vector_time,ana4_epochs(i,:),'k');hold on;
    % %             plot(vector_time,co2(i,:),'r:');
    % %             legend('TDT-pressure', 'SL-pressure');
    % %             title(['epoch #' num2str(i)]);
    % %             pause(2)
    % %         end
    % %     end

    % %     %%ana5
    % %     ana5 = ecog.smartlab.tdt.o2;
    % %     [ana5_epochs] = BU_epochSignal(ana5, fs1, TDTtimestamps, pre, pos);
    % %     o2ori = squeeze(SL_epochsplus(3,:,:));
    % %     o2 = zscore(o2ori);
    % %     ana5 = zscore(ana5_epochs);
    % %
    % %     if plotting_tr
    % %         for i = 1:76;clf;
    % %             plot(vector_time,ana5_epochs(i,:),'k');hold on;
    % %             plot(vector_time,o2(i,:),'r:');
    % %             legend('TDT-pressure', 'SL-pressure');
    % %             title(['epoch #' num2str(i)]);
    % %             pause(2)
    % %         end
    % %     end


    %% Summary plot of epochs - Gold Standard!!
    %  respiration signal recorded in TDT & Pressure signal recorded in SL
    plotting_tr = 1;
    ntrs = numel(timestampsSL);
    if plotting_tr && isfield(ecog.epochs,'resp')
        figure
        for i = 1:ntrs;clf;
            x = ecog.epochs.time{i,1};
            y1 = zscore(ecog.epochs.resp{i}');
            y2 = zscore(ecog.epochs.pressure{i});
            %y3 = zscore(ecog.epochs.flow{i});

            plot(x,y1,'k');hold on; % resp recorded on tdt (belt)
            plot(x,y2,'r');hold on;

            title(['epoch #' num2str(i)]);
            pause; clf
        end
    end

    % %% Sam added summary plots for the three additional analog signals fo vccp protocol
    % if plotting_tr
    %     figure
    %     for i = 1:ntrs; clf;
    %         subplot(411)
    %         plot(ecog.epochs.time, demean(ecog.epochs.pressure{i}),'r'); hold on; % FROM SAMRTLAB DIRECTLY
    %         plot(ecog.epochs.time, demean(ecog.epochs.pressure_tdt{i}),'b:'); % FROM TDT
    %         legend('Pressure (TDT sl)', 'Raw Pressure (TDT)');
    % 
    %         subplot(412)
    %         plot(ecog.epochs.time, demean(ecog.epochs.co2{i}),'r'); hold on; % FROM SAMRTLAB DIRECTLY
    %         plot(ecog.epochs.time, demean(ecog.epochs.co2_tdt{i}),'b:'); % FROM TDT
    %         legend('CO2 (TDT sl)', 'Pct CO2 (TDT)');
    % 
    %         subplot(413)
    %         plot(ecog.epochs.time, demean(ecog.epochs.flow{i}),'r'); hold on;
    %         plot(ecog.epochs.time, demean(ecog.epochs.gsm_mixer{i}),'b:');
    %         legend('Flow (TDT sl)', 'Outflow Rate (TDT)');
    % 
    %         title(['epoch #' num2str(i)]);
    %         pause;
    %     end
    % end

    %% Sam added above block for additional analog signals foR vccp protocol - end


    %% summary plot - other signals

    if plotting_tr
        figure
        for i = 1:ntrs;clf;

            % subplot(411)
            % plot(ecog.epochs.time,ecog.epochs.pressure{i},'r');hold on;
            % plot(ecog.epochs.time,ecog.epochs.pressure_tdt{i}*4,'k');hold on;
            % if isfield(ecog.epochs,'inhonsets_tdt') && ~isempty(ecog.epochs.inhonsets_tdt)
            %     vline(ecog.epochs.inhonsets_tdt{1},'k:');
            % end
            % if isfield(ecog.epochs,'inhonsets') && ~isempty(ecog.epochs.inhonsets)
            %     vline(ecog.epochs.inhonsets{1},'r:');
            % end
            % legend('SL-pressure', 'SL-TDT-pressure');

            subplot(412)
            plot(ecog.epochs.time{i,1},demean(ecog.epochs.pressure{i}),'r');hold on;
            plot(ecog.epochs.time{i,1},demean(ecog.epochs.co2{i}),'r:');hold on;
            legend('Pressure (SL)', 'CO2(SL)');
            vline(ecog.epochs.inhonsets{1},'r');

            % subplot(413)
            % plot(ecog.epochs.time,ecog.epochs.co2{i},'r');hold on;
            % plot(ecog.epochs.time,ecog.epochs.co2_tdt{i}*4,'k');hold on;
            % legend('CO2 (SL)', 'CO2(SL-TDT)');
            % vline(ecog.epochs.inhonsets{1},'r');


            subplot(414)
            plot(ecog.epochs.time{i,1},demean(ecog.epochs.pressure{i}),'r');hold on;
            plot(ecog.epochs.time{i,1},demean(ecog.epochs.o2{i}),'r:');hold on;
            legend('Pressure (SL)', 'O2(SL)');

            % figure
            % subplot(411)
            % plot(ecog.epochs.time,demean(ecog.epochs.pressure{i}*10),'r');hold on;
            % plot(ecog.epochs.time,demean(ecog.epochs.flow{i}),'r:');hold on;
            % legend('Pressure (SL)', 'Flow(SL)');

            %% 
            % subplot(412)
            % plot(ecog.epochs.time,demean(ecog.epochs.pressure{i}*1000),'r');hold on;
            % plot(ecog.epochs.time,demean(ecog.epochs.volume{i}),'r:');hold on;
            % legend('Pressure (SL)', 'Volume(SL)');
 
            %plot(ecog.epochs.time,ecog.epochs.co2{i},'r:');
            %vline(ecog.epochs.inhonsets_tdt{1},'k');
            %plot(ecog.epochs.time,ecog.epochs.inhdur{i},'k');

            title(['epoch #' num2str(i)]);
            close
            pause;%(2)
        end
    end

   %save([params.directoryOUT params.filename '_ecog_epoch_smartlab_o.mat'], 'ecog', '-v7.3');
   save([params.directoryOUT params.filename '_ecog_epoch_smartlab.mat'], 'ecog', '-v7.3');

end

%% STEP 4.5
% ADD processed SL DATA (SLtable.mat) to 'ecog' structure.
% Align signals to SINK pulses - the common digital input pulse sent to both amps
% Sink pulses are usually 'ringExp' (in RRST/RRET) or 'breathswitch' (balonX) in P100
% Epoch all signals (ex. time window -15:45) 
% Save ready2rock data 'ecog' variable (ex. B48_P100_ecog_epoch_smartlab.m
% COPY VARIABLE TO M2-CHIP AND RENAME IT ('NS195_RRET.M', 'NS196_B48_P100.M')

% JA added: preprocess ET data [should go after Step 4 (after ecog is epoched)]
% [1] make sure et setup is correct in params
% [2] fill out BXX_sink_et in log file
% make this new sheet, delete unused TTL rows, and fill in TTLs (3:end-1; these are first 2, then
% ring expansions)
% [3] add et to ecog

addETnow = 1;
if addETnow

% load([params.directoryOUT params.filename '_ecog_epoch_smartlab.mat']) % loads trial-epoched TDT, SL data

% manually add ET timestamps to sink spreadsheet
trialI = 7; % manually add et TTLs
for i = 1 %trialI = 1:length(params.et.currentBlockRecs)
    recNum = params.et.currentBlockRecs(trialI);
    trNum = params.et.currentBlockTrs(trialI);
    etRecDir = [params.et.dir 'rec' num2str(recNum) '_tr' num2str(trNum)];

    eventTable = readtable([etRecDir '/events.csv'], 'Delimiter',',');
    eye_stateTable = readtable([etRecDir '/eye_state.csv'], 'Delimiter',',');

    et_T0 = eye_stateTable.ts;
    et_event_ts = eventTable.timestamp_ns_/1000000000; % *** subtract out t1 ?
    %et_TTL_ts = et_event_ts(3:end-1);
    et_TTL_ts = et_event_ts(1:2);
    buddy = diff(et_TTL_ts); % manually inspect this
end

%% 
addETdata = 1;
if addETdata
if isfield(params.et,'sinktimes') && ~isempty(params.smartlab.sinktimes)
    % load sink timestamps from excel file (in secs)
    idx = 1:strfind(params.et.sinktimes, '.xlsx');
    filenamenow = [params.et.sinktimes(idx) 'xlsx']; % excell fullname including sheet
    %sheet = params.et.sinktimes(idx(end)+5:end); %initially, included
    %'/'before sheet name gave error
    sheet = params.et.sinktimes(idx(end)+6:end);
    SinkTable = readtable(filenamenow, 'Sheet', sheet); % load excell sheet with sink timestamps (in secs): 'SL_time', 'TDT_time'
    % provide titles if not picked up
    % SinkTable.Properties.VariableNames(1) = "SL_time";
    % SinkTable.Properties.VariableNames(2) = "TDT_time";
    % SinkTable.Properties.VariableNames(3) = "ET_time";
    % SinkTable.Properties.VariableNames(4) = "IDs";
    % SinkTable.Properties.VariableNames(5) = "Trial";
    % SinkTable.Properties.VariableNames(6) = "Pulse_num";
    % extract timestamps
    Cnames = SinkTable.Properties.VariableNames'; % 'SL_time', 'TDT_time'
    sinkpulsesTDT = SinkTable.TDT_time;
    sinkpulsesET  = SinkTable.ET_time;
end
if strcmp(params.task,'VCCP')
    pre = 60; pos = 160; % select time window around event - 'ringexp1' is time 0s ***JA: !should be same as in Step 4!
    numringexp = 20; % ring expands 10 times [***JA: different # of expansions on different trials...]
    sinkID = 'trialstart';
    stopID = 'trialstop';
end

% original code only this, sam added the following block to reflect the use
% of new function used to cut trials based on the TDT trials instead of pre
% and post seconds 
% index = find(strcmp(SinkTable.IDs, sinkID)==1); % select event type (IDs)
% timestampstdt = sinkpulsesTDT(index);
% timesptampsET  = sinkpulsesET(index);

index = find(strcmp(SinkTable.IDs, sinkID)==1); % select event type (IDs)
stop_index = find(strcmp(SinkTable.IDs, stopID)==1);
timestampstdt = sinkpulsesTDT(index);
timestampsET  = sinkpulsesET(index);
timestampstdt_stop = sinkpulsesTDT(stop_index);
timestampsET_stop = sinkpulsesET(stop_index);

TDT_trial_dur = timestampstdt_stop - timestampstdt;

ntrs = length(ecog.epochs.eeg);
ecog.epochs.et = [];

% sam added : when rec_tr ET are not recorded for all TDT trials for the
% block
% Make sure mapping vectors match number of EEG trials
% Fill missing ET trial(s) with NaN
params.et.currentBlockRecs = [7 8 9 10 11 12 13];
params.et.currentBlockTrs  = [8 9 10 11 12 13 14];

% ending et mapping 

for trialI = 1:ntrs

    recNum = params.et.currentBlockRecs(trialI);
    trNum = params.et.currentBlockTrs(trialI);

    missingTrial = isnan(recNum) || isnan(trNum) || isnan(timestampsET(trialI)) || isnan(timestampsET_stop(trialI));
    if missingTrial
        continue
    end

    etRecDir = [params.et.dir 'rec' num2str(recNum) '_tr' num2str(trNum)];
    blinksTable = readtable([etRecDir '/blinks.csv']);
    eventTable = readtable([etRecDir '/events.csv'], 'Delimiter',',');
    eye_stateTable = readtable([etRecDir '/eye_state.csv'], 'Delimiter',',');
    if log10(eye_stateTable.ts(1))>18 % in ns
        eye_stateTable.ts = eye_stateTable.ts/1000000000; % convert to ms
    end

    isBlink = zeros(size(eye_stateTable, 1), 1);
    for blinkI = 1:size(blinksTable, 1)
        tI1 = findnearest(blinksTable.startTimestamp_ns_(blinkI)/1000000000, eye_stateTable.ts);
        tI2 = findnearest(blinksTable.endTimestamp_ns_(blinkI)/1000000000, eye_stateTable.ts);
        isBlink(tI1:tI2) = 1;
    end
    eye_stateTable.isBlink = isBlink;

    % trial start relative to first ET sample
    timestamp1 = timestampsET(trialI) - eye_stateTable.ts(1);
    timestamp1_stop = timestampsET_stop(trialI) - eye_stateTable.ts(1);

    column_names = eye_stateTable.Properties.VariableNames';% columns titles (Time,Flow,Volumen,Pressure,...)
    for j = 1:numel(column_names)
        column_name = column_names{j};
        column_values = eye_stateTable.(column_name);
        
        %% upsample freq to 500hz and epoch
        oldfs = params.et.fs; newfs = ecog.epochs.fs;
        if isnan(timestampsET(trialI)) || isnan(timestampsET_stop(trialI))
            ecog.epochs.et.(column_name){trialI,1} = [];
        else
            % trial exists - call BU_EpochDataTrials
            [epochsnow] = BU_EpochDataTrials(column_values, oldfs, newfs,timestamp1, timestamp1_stop, TDT_trial_dur(trialI));
            ecog.epochs.et.(column_name){trialI, 1} = epochsnow{1};
        end
    end
end

% plot resulting epoched ET data
T = ecog.epochs.time;
epochFields = fieldnames(ecog.epochs.et);
for trialI = 1:ntrs
    nexttile
    isEpochedVar = zeros(size(epochFields));
    for fieldI = 2:length(epochFields)
        buddy = ecog.epochs.et.(epochFields{fieldI});
        if iscell(buddy) && (min(size(buddy)==[ntrs 1])==1) && length(buddy{trialI}) == length(T) && ~strcmp(epochFields{fieldI}, 'eeg')
            plot(T, buddy{trialI})
            hold on
            isEpochedVar(fieldI) = 1;
        end
    end
    legend(epochFields(find(isEpochedVar)))
end

save([params.directoryOUT params.filename '_ecog_epoch_smartlab_et.mat'],'ecog', '-v7.3');
end
%% 

% JA added: instead of step 5 (manually fine-tune inh ons/offs) use this
% first, make sure that trialsinfo sheet is all filled out
% and, make sure to add "-<area>" to elec correspondence sheet
load([params.directoryOUT params.filename '_ecog_epoch_smartlab.mat'],'ecog'); %LOAD FINAL ET FILE IN NORMAL CASES, for NS215, NO ET DATA
epochs = ecog.epochs;
chansinfo = readtable(ecog.params.labelfile, 'Sheet','TDT');
params = ecog.params;
trialsinfo = readtable(params.smartlab.sinktimes(1:end-13), 'Sheet', [params.filename(1:end) '_trialsinfo']);
% JA *** !uncomment below if the file is not named correctly!
%directoryclean.smartlabfile = [directoryclean.smartlabfile(1:end-8) params.filename(1:3) directoryclean.smartlabfile(end-8:end)]; % JA added, otherwise no distinction between blocks
save(directoryclean.smartlabfile, 'epochs', 'trialsinfo', 'chansinfo', 'params', 'directoryclean')

end

%% --------------------------------------------------------------------------------------------
if Find_ChanswithRRBOs

    % RRBOs analyses two steps 
    % 1) epoching all data types around inh onsets (RRBOs) 2) computing stats
    % include ALL inhs (loaded and non-loaded)
    % epoch time-window: 2s before inh onset and 3s after
    % include inhs from ALL trials (bad ones also as we need as many inhs as possible for RRBO analyses)
    % input data should be 'epochs' structure with TDT/SL fields & inspected inh onset times

    clear Step* Find_*
    %directoryclean
    %filename
    disp(['LOADING...' directoryclean.smartlabfile])

    for f = 1:length(filename)

        %% LOAD 'EPOCHS'
        % JA *** !uncomment below if the file is not named correctly!
        %directoryclean.smartlabfile = [directoryclean.smartlabfile(1:end-8) params.filename(1:3) directoryclean.smartlabfile(end-8:end)]; % JA added, otherwise no distinction between blocks

        % 'epochs' structure with TDT/Natus + SL fields AND inspected inh onset times
        load(directoryclean.smartlabfile); % load struct
        
        %clearvars -except epochs params chansinfo trialsinfo filename directoryclean
      
        %% ALIGN EPOCHS TO 'INHONSETS' 
        % epochs were originally alligned to 'sinkID' 
        % align epochs to 'InhOnset'--> the first inh onset after load/occlusion

        if strcmp(params.task,'P100') % align only in P100 because time 0 should be inh onset (for RRST/RRET time 0 should be ringexp#1)
            if ~strcmpi(epochs.sinkID,'InhOnset')  % check if epochs aligned already
                delayoffset = 0;%-0.2;
                plotting = 0;
                trialsinfo; %[]empty if not provided 

                [epochs,trialsinfo] = BU_realign_epochs (epochs,delayoffset,trialsinfo,plotting);
            end
        end
       
        %% RRBO ANALYSES - STEP 1
        %  1- filter, bipolarize, detect IEDs
        %  2- add fields to epochs: 'epochs.IEDs', 'epochs.eegbip',
        %  3- ask user to resave new 'epochs' to file (NEEDS MANUAL SAVING)

        %% filter, bipolarize, detect IEDs

        % message = 'FILTER, BIPOLARIZE, DETECT IEDS ?';
        % user_response = prompt_save_dialog(message);
        user_response = 1; % JA overwrite

        if user_response
            if isfield (epochs,'eegbip_IEDs') && isfield(epochs,'eegbip') % check if bipolization done; check if IED done
                disp('BIPOLARIZATION AND IEDs ALREADY DONE');pause
            end
            if isfield (epochs,'eeg_IEDs') && isfield(epochs,'eeg') %  check if IED done
                disp('IEDs ALREADY DONE');pause
            end

            % select params for epochs:
            paramsepochs = [];
            paramsepochs.bipolarize = 0; %[1]bipolorize to near WM; [0]monopolar
            paramsepochs.detectIEDs_plus = 1;
            paramsepochs.filter_bandpass = [];%[0.1 95]; %[]leave empty if not neeeded  %some SL data is noisy at higher freq
            paramsepochs.plotting = 0; % plot IEDs
            %paramsepochs.reference = 'none';
            
           % add electrodes info
            paramsepochs.chansinfo = chansinfo;
            
            % make destination folder (to save results)
            savefigfolder = {};
            % filenameshort = strtok(filename{1}, '_');
            filenameshort = strtok(filename, '_'); %samruddhi added, indexing error with prev line

            filePath = directoryclean.smartlabfile;
            [folderPath, ~, ~] = fileparts(filePath);% Extract the directory part of the file path
            %newDirPath = fullfile(folderPath, 'fig', 'beh'); % Define the new directory path by appending 'fig/beh' to the existing path
            
            if paramsepochs.bipolarize
                %savefigfolder = ['/Users/jherreroru/Documents/data/Dyspnea/RRST/' params.subjectID '/report/figs/RRBO_bip/']; %[]do not save
                %savefigfolder = ['/Users/jherreroru/Documents/data/Dyspnea/P100/' params.subjectID '/' filenameshort '/RRBO_bip/']; %[]do not save
                savefigfolder = fullfile(folderPath, 'RRBO_bip'); % Define the new directory path by appending 'RRBO_bip' to the existing path
            else
                %savefigfolder = ['/Users/jherreroru/Documents/data/Dyspnea/P100/' params.subjectID '/' filenameshort '/RRBO/']; %[]do not save
                savefigfolder = fullfile(folderPath, 'RRBO');
            end
            
            if ~exist(savefigfolder, 'dir')
                mkdir(savefigfolder);
            end
            
            paramsepochs.savefigfolder = savefigfolder;
            
            clear newDirPath filePath folderPath

            BU_RRBO_addfields (epochs,paramsepochs);% DO IT! % JA: have to go into this, manually save at end
            disp('step 1 done')
        end

        %% RRBO ANALYSES - STEP 2
        % epoch all Inhs into 'epochs.RRBO'
        % input data are 'epochs'--> trials with multiple inhs each
        % remove Inhs with spks (if 'discardIEDs' is not empty)
        % remove Inhs with loads (if 'discardLoads' equals 1)
        % add field 'trialsinfo.inhsinfo' with info about each inh (ex,load vs no-load)

        %load(directoryclean.smartlabfile) % JA added
        load('/home/david/Documents/RRET/patients/NS217/vccp/NS217_VCCP_B30_P1.mat');
        %sam commented above line, check later if needed to uncomment
        %%sam commented for NS215

        % SELECT BIPOLAR OR MONOPOLAR
         disp(['BIPOLARIZING NOW =' num2str(epochs.paramsepochs.bipolarize)]);
         pause

        if  isfield(epochs,'paramsepochs') 
            if epochs.paramsepochs.bipolarize == 1; paramsepochs.signalref = 'eegbip'; end
            if epochs.paramsepochs.bipolarize == 0; paramsepochs.signalref = 'eeg'; end
        end

        %paramsepochs = [];
        paramsepochs.bastime       =  -8; %used to be -6 % desired baseline time (pre InhOnset); use at least -4 in RRST/RRET
        paramsepochs.postime       =  8; % used to be 6 % desired post time (post InhOnset); use at least -4 in RRST/RRET
        paramsepochs.extratime     =  12; % last time point allowed to include last inh after ring #4 (ex. 4s )
        paramsepochs.discardIEDs   =  0; % [1] discard epochs with IEDs around inh onset (-1:1); [0] do not discard -include all inhs even those with IEDs - add an IED matrix with index inhs with no spks with 1s
        paramsepochs.Loads         =  2; %[0,1,2] % [2]include ALL inhs (non-loaded & loaded);[0]only inhs WITHOUT load/occl (non-loaded); [1]only inhs WITH loads/occlusions;
        paramsepochs.filtfilt      = [];%[.01 150] filter eeg a bit more (avoid pump noises);[]do not filter

        paramsepochs.savefigfolder = [];%[epochs.paramsepochs.savefigfolder '/fig']; % [] do not save
        paramsepochs.savefolder    = [paramsepochs.savefigfolder '/'];%[epochs.paramsepochs.savefigfolder '/'];% '/TF_bip/']; % [] do not save
        paramsepochs.task          = 'P100';%'P100';

        if exist('trialsinfo','var') % info about trials (eg, for RRST experiment)
            paramsepochs.trialsinfo = trialsinfo;
        end
        % paramsepochs.signalref = 'eeg'; % JA added, but it should already
        % exist (if you ran IEDs)

        [RRBO] = BU_RRBO_epoch_clean (epochs, paramsepochs); % DO IT!
       
        % save ([params.directoryOUT 'RRBO'], 'RRBO','chansinfo','trialsinfo'); %name variable (e.g., NS195_RRST_RRBO.mat)
        save([savefigfolder '/' filename{1}(1:3) '_RRBO.mat'], 'RRBO','chansinfo', 'trialsinfo', 'params', '-v7.3'); % JA replaced above with this

       %% ******* YOU CAN STOP HERE AND JUMP TO 'ALLPATS_RRST_RESULTS.m '  *******

        %% ----------------------------------------------------------------------------------------------------------------------------
        %% ----------------------------------------------------------------------------------------------------------------------------

        %% RRBOs analyses

        %% load selected electrodes from params file

        if isfield(params,'rrbo_siglabels') && ~isempty(params.rrbo_siglabels)
            RRBO.rrbo_siglabels = params.rrbo_siglabels;
        end

        if isfield(params,'rrst_siglabels') && ~isempty(params.rrst_siglabels)
            RRBO.rrst_siglabels = params.rrst_siglabels;
        end

        if isfield(params,'rret_siglabels') && ~isempty(params.rret_siglabels)
            RRBO.rret_siglabels = params.rret_siglabels;
        end

        if isfield(params,'p100_siglabels') && ~isempty(params.p100_siglabels)
            RRBO.p100_siglabels = params.p100_siglabels;
        end

        %% select specific electrodes (& formatting options)

        selidx = 1:numel(RRBO.eeg); % all chans selected

        % Define the options
        options = {'Include ALL electrodes', 'Include *** electrodes'};
        % Display the GUI dialog
        [selectedIndex, tf] = listdlg('PromptString', 'Choose an option:', ...
            'SelectionMode', 'single', ...
            'ListString', options);

        if selectedIndex == 2

            % Define field options
            fieldOptions = {'rrbo_siglabels', 'rrst_siglabels', 'rret_siglabels', 'p100_siglabels'};

            % Display GUI
            [selectedIndex, tf] = listdlg('PromptString', 'Choose electrodes:', ...
                'SelectionMode', 'single', ...
                'ListString', fieldOptions, ...
                'Name', 'Field Selection');

            % Check if user made a choice
            if tf % if user made a choice
                selectedField = fieldOptions{selectedIndex};

                % Check if the selected field is empty
                if isempty(RRBO.(selectedField))
                    %warndlg('The selected field is empty.', 'Empty Field');
                else
                    selchans = RRBO.(selectedField);
                end
            else % if the user closed the dialog or pressed "Cancel"
                %warndlg('User did not make a selection.', 'No Selection Made');
                return;
            end

            % Format levels
            if isfield(epochs,'labels_bip') && ~isempty(epochs.labels_bip)
                % only take elec#1 (from the pair)
                monopolarLabels = {};
                for chan = 1:length(epochs.labels_bip)
                    idx = findstr(epochs.labels_bip{chan},'-');
                    monopolarLabels{chan,1} = epochs.labels_bip{chan}(1:idx-1);
                end
                % select sigerps chans
                selidx = find(ismember(monopolarLabels,selchans)==1);
            end
        end

        %% APPLY ELECTRODE SELECTION

        if isfield(RRBO,'eeg') && ~isempty(RRBO.eeg)
            eeg = RRBO.eeg(selidx);
        end

        %% COMPUTE DOMINANT FRQ FOR EVERY ELECTRODE

        RRBO.domfrq = [];
        % Define field options
        fieldOptions = {'compute dominant frq', 'close fig'};

        % Display GUI
        [selectedIndex, tf] = listdlg('PromptString', 'Choose electrodes:', ...
            'SelectionMode', 'single', ...
            'ListString', fieldOptions, ...
            'Name', 'Field Selection');

        if tf % if user made a choice
            selectedField = fieldOptions{selectedIndex};

            % Check if the selected field is empty
            if strcmp(selectedField,'compute dominant frq')

                RRBO.domfrq = [];

                % declare params
                twin = [-0.5 1]; % ** define time-win for power peak **
                plotting = 0;

                eeg;
                time_vector = RRBO.time;
                savefolder = paramsepochs.savefigfolder;
                result = gpt_find_dominant_frq2 (eeg, time_vector, twin, savefolder,plotting);

                RRBO.domfrq = result;% collect
                clear result

            end
        end

        %% fomat data for next step 
        % format data and delete repeated data
        % take only first cell because all are equal 

        if isfield(RRBO,'pressure') && ~isempty(RRBO.pressure)
            pressure = RRBO.pressure(selidx);
            if paramsepochs.discardIEDs == 0
                pressure = pressure(1); % {samples * sweeps} ... delete repeated data
            end
        end

        if isfield(RRBO,'flow') && ~isempty(RRBO.flow)
            flow = RRBO.flow(selidx);
            if paramsepochs.discardIEDs == 0
                flow = flow (1); % {samples * sweeps}
            end
        end

        if isfield(RRBO,'volume') && ~isempty(RRBO.volume)
            volume = RRBO.volume(selidx);
            if paramsepochs.discardIEDs == 0
                volume = volume (1); % {samples * sweeps}
            end
        end

        if isfield(RRBO,'emg') && ~isempty(RRBO.emg)
            emg = RRBO.emg(selidx);
            if paramsepochs.discardIEDs == 0
                emg = emg (1); % {chans * samples * sweeps}
            end
        end

        if isfield(RRBO,'resp') && ~isempty(RRBO.resp)
            resp = RRBO.resp(selidx);
            if paramsepochs.discardIEDs == 0
                resp = resp(1); % {chans * samples * sweeps} (or samples sweeps if only 1 chan)
            end
        end

        if isfield(RRBO,'co2') && ~isempty(RRBO.co2)
            co2 = RRBO.co2(selidx);
            if paramsepochs.discardIEDs == 0
                co2 = co2(1); % {samples * sweeps}
            end
        end

        if isfield(RRBO,'o2') && ~isempty(RRBO.o2)
            o2 = RRBO.o2(selidx);
            if paramsepochs.discardIEDs == 0
                o2 = o2(1); % {samples * sweeps}
            end
        end

        if isfield(RRBO,'cranial') && ~isempty(RRBO.cranial)
            cranial = RRBO.cranial(selidx);
            if paramsepochs.discardIEDs == 0
                cranial = cranial(1);% {chans * samples * sweeps}
            end
        end

        time             = RRBO.time;

        totalinhs        = RRBO.totalinhs;

        labels_bip       = epochs.labels_bip(selidx);
        labels_bip_FSurf = epochs.labels_bip_FSurf(selidx);
        
        %% check trialinfo

        if exist('trialsinfo','var')
            inhsinfo = RRBO.inhsinfo;

            if isfield(inhsinfo,'trnum') && ~isempty(inhsinfo.trnum)
                inhsinfo.trnum = inhsinfo.trnum(selidx);
            end

            if isfield(inhsinfo,'load') && ~isempty(inhsinfo.load)
                inhsinfo.load = inhsinfo.load(selidx);
            end

            if isfield(inhsinfo,'loaddur') && ~isempty(inhsinfo.loaddur)
                inhsinfo.loaddur = inhsinfo.loaddur(selidx);
            end

            if ~isempty(inhsinfo.criteria)
                inhsinfo.criteria   = inhsinfo.criteria(selidx);
            end

            if ~isempty(inhsinfo.good)
                inhsinfo.good       = inhsinfo.good(selidx);
            end

            if ~isempty(inhsinfo.intensity)
                inhsinfo.intensity  = inhsinfo.intensity(selidx);
            end

            if ~isempty(inhsinfo.confidence)
                inhsinfo.confidence = inhsinfo.confidence(selidx);
            end

            if ~isempty(inhsinfo.breathless)
                inhsinfo.breathless = inhsinfo.breathless(selidx);
            end

            if ~isempty(inhsinfo.keypress)
                inhsinfo.keypress = inhsinfo.keypress(selidx);
            end

            if isfield(inhsinfo,'ringnum') && ~isempty(inhsinfo.ringnum)
                inhsinfo.ringnum = inhsinfo.ringnum(selidx);
            end

            %% clean-up inhsinfo (clear empty fields)

            new_inhsinfo = struct(); % Create an empty structure for the new data

            fields = fieldnames(inhsinfo); % Get the names of all fields in the original structure

            for i = 1:length(fields)
                field_name = fields{i};
                if ~isempty(inhsinfo.(field_name))
                    new_inhsinfo.(field_name) = inhsinfo.(field_name); % Copy the field if it's not empty
                end
            end
            inhsinfo = new_inhsinfo;
            inhsinfo.totalinhs = totalinhs;
        end

        %% prepare input vv
        dataIN = [];
        if exist('eeg','var');dataIN.eeg = eeg;end

        if exist('eeg_IEDs','var');dataIN.eeg_IEDs = eeg_IEDs;end

        if exist('cranial','var');dataIN.cranial = cranial; end

        if exist('emg','var');dataIN.emg = emg; end

        if exist('resp','var');dataIN.resp = resp; end

        if exist('time','var');dataIN.time = time;end

        if exist('labels_bip','var');dataIN.labels_bip = labels_bip;end

        if exist('labels_bip_FSurf','var');dataIN.labels_bip_FSurf = labels_bip_FSurf;end

        if exist('pressure','var');dataIN.pressure = pressure;end

        if exist('o2','var'); dataIN.o2 = o2; end

        if exist('co2','var'); dataIN.co2 = co2; end

        if exist('flow','var'); dataIN.flow = flow; end

        if exist('volume','var'); dataIN.volume = volume; end


        %% ******* COMPUTE ERPs (Inh-locked) ******* 
        % after visual inspection, user needs to write selected channels to
        % epochs.rrbo_sigerpslabels and params.rrbo_siglabels

        % select params
        paramsnow = [];
        paramsnow.basecorr = [-1.5 -.9]; %  [-1 -.5] baseline correction; [] no correction (will not compute stats) ; 0ms is pressure onset
        paramsnow.statswin = [-.6 1.5]; % window for overall significance (one p-val for the whole time-window)
        paramsnow.xlimshow = [-1.5 3]; % xlim for fig
        paramsnow.savefigfolder =  paramsepochs.savefigfolder;

        % do it
        sig_erp = BU_RRBO_erps (dataIN,paramsnow,inhsinfo); % all inhs together

        %% (NEEDS WORK)
        % compute RRBO - FrqBands
        % TF analysis specfic freq bands using hilbert, after vis inpspection write selected channels to epochs.rrbo_sigfrqbandlabels

        % [RRBO.sig_frqband] = BU_RRBO_frqbands (eeg,pressure,time,labels_bip,labels_bip_FSurf,...,
        %                paramsepochs.savefigfolder,paramsepochs.basecorr,totalinhs,inhsinfo);
       
        %% ---------------------------------------------------------------------------
        %  *******  COMPUTE RRBO - TF SPECTRUM ******* 
        % compute power/itc wavelets for all chans or selected channels

        % select params
        paramsnow = [];
        paramsnow.timeframebas = []; %[-1 -.5] baseline correction; [] no correction (0ms is inh onset based on pressure signal)
        paramsnow.timeframe    = [dataIN.time(1)+.1 dataIN.time(end)-.1];%[-1.5 3];%[-1 2];%[dataIN.time(1)+.1 dataIN.time(end)-.1];% smaller that original window size to avoid edge effects
        paramsnow.frqlo        = .01;
        paramsnow.frqhi        = 150;
        paramsnow.paramsori    = params; 

        % Check if folder exists (create it otherwise)
        savefolder = [paramsepochs.savefigfolder 'TF_allchans1/'];
        if ~exist(savefolder, 'dir')
            mkdir(savefolder);
        end
        paramsnow.savefigfolder = savefolder;

        %% DO IT
        %sig_tf = BU_RRBO_tfpsd_beta (dataIN,paramsnow,inhsinfo);%modify this to save single chans -add missing fields (resp,flow,emg
        sig_tf = BU_RRBO_tfpsd (dataIN,paramsnow,inhsinfo);


        %% same as above but saves eaach field to a separate file
        % ---> also adds all other fields (resp,flow,emg,..)
        sig_tf_plus = BU_RRBO_tfpsd_addfields (dataIN,paramsnow,inhsinfo);


        %% COMPUTE RRBO-FrqBand
        % extracts specific fbands from the above analyses (wavelets)
        % averages itc/power for each frqband
        % create new variable
        % ---> needs output from 'sig_tf' above

        % select params
        sig_tf.timeframebas = [];%[-1.2 -0.8];%[];
        sig_tf.timeframe = [-2 2];
        plotting = 1;

        % prepare
        clearvars -except sig_tf paramsepochs inhsinfo plotting
        savefolder = [paramsepochs.savefolder];
        if ~exist(savefolder, 'dir')% Check if the folder exists; if not, create it
            mkdir(savefolder);
        end
        load([savefolder 'sig_domfrq']);

        % do it
        sig_tfb = BU_avg_freq_bands(sig_tf,sig_domfrq,savefolder,plotting);
        %----> save sig_tfb to file and continue (ex. 'sig_tfb_power.mat')
        % if vv already saved, start debugging from here
        % if exists('sig_tf')
        %     save ....
        %   else
        %     load ...
        % end

        %% LOAD vs noLOAD - FrqBand
        % sig_tf.m is huge!
        % - extract only necessary fields
        % - select only ** electrodes

        clearvars -except sig_tfb paramsepochs


        plotting = 1;
        eegfb = sig_tfb.pow_st_fb;%{chans}(trs*fbs*samps)
        pressure = sig_tfb.pressure_st;%{chans}(trs*samps)
        inhsinfo = sig_tfb.inhsinfo;

        % check that the 2 vv are equal (meaning that they were extracted properly
        pressure1 = pressure{1};
        savefolder = [paramsepochs.savefolder 'selchan/singlevv/'];
        load([savefolder 'pressure_st']);
        pressure2 = pressure_st{1};

        %% respiration or flow
        load([savefolder 'resp_st']);


        %% params (select)
        paramsnow = [];
        paramsnow.comparewhat = 'load'; % 'loadcorrect';
        paramsnow.timeframebas = [];%paramsepochs.basecorr;
        paramsnow.plotting = 1;
        paramsnow.task = 'P100';

        %% additional params
        savefolder = [paramsepochs.savefolder 'fig_p100_2/'];
        if ~exist(savefolder, 'dir') % Check if the folder exists; if not, create it
            mkdir(savefolder);
        end
        paramsnow.savefigfolder = savefolder; clear savefolder
        paramsnow.time = single(sig_tfb.time);
        paramsnow.labels_bip = sig_tfb.labels;
        paramsnow.labels_bip_FSurf = {};%sigtf.labels_bip_FSurf;
        paramsnow.totalinhs = max(cellfun(@(x) size(x, 1), inhsinfo.trnum));
        paramsnow.frqbands = sig_tfb.fblim;

        [RRBO.rrst_sigerpfb] = BU_RRST_erps_fb (eegfb,pressure,resp_st,inhsinfo,paramsnow)

        %% compute ERPs: LOAD vs. NO-LOAD (RRST/RRET/P100 erp (all or selected channels)

        comparewhat = 'load'; %'load'; 'load-correct';

        [RRBO.rrst_sigerp] = BU_RRST_erps (eeg,pressure,flow,time,labels_bip,labels_bip_FSurf,...,
            paramsepochs.savefigfolder,timeframebas,paramsepochs.task,...,
            totalinhs,inhsinfo,comparewhat);

        %% compute ERPs freq-band specific anal : LOAD vs. NO-LOAD (RRST/RRET/P100 erp (all or selected channels)
        % [RRBO.rrst_sigerp] = BU_RRST_erps_fband (eeg,pressure,flow,time,labels_bip,labels_bip_FSurf,...,
        %     paramsepochs.savefigfolder,timeframebas,paramsepochs.task,...,
        %     totalinhs,inhsinfo,comparewhat);


        %% compute TF: LOAD vs. NO-LOAD (RRST/RRET/P100) (all or selected channels)

        timeframe = [-1.5 1.5];
        timeframebas = paramsepochs.basecorr;
        frqlo = 0.05;
        frqhi = 60;
        %% CREATE THIs NEW PROGRAM - modify BU_RRBO_tfpsd to compare 2 conditions
        [RRBO.p100_sigtf] = BU_RRST_tfpsd (eeg,pressure,time,labels_bip,labels_bip_FSurf,...,
            paramsepochs.savefigfolder,timeframebas,timeframe,frqlo,frqhi,...,
            totalinhs,inhsinfo,comparewhat);


    end
end

%% when co2 in tdt wasn't in analog but in a different stream - Resp in the tdt block
% After step3 run this
% filename = char("/home/david/Documents/RRET/patients/NS217/vccp/tdt/VCCP_NS217_B30_P1");
% tdt = TDTbin2mat(filename);
% 
% %load("/home/david/Documents/RRET/patients/NS217/vccp/mat/VCCP_NS217_B30_P1_ecog_epoch.mat");
% if isfield (params.smartlab.tdt,'co2')  && ~isempty(params.smartlab.tdt.co2)
% 
%     Sdelay  = params.smartlab.tdt.co2(2);
%     Stype   = 'analog';
%     Schan   = params.smartlab.tdt.co2(1);% analog chan were CO2 was recorded
%     Ssignal = tdt.streams.Resp.data; %sam changed
%     Sfs     = tdt.streams.Resp.fs;   %sam changed
%     Slabel   = 'pct CO2 smartlab';
%     correctdelay = 1;
%     if Sdelay ~= 0 && correctdelay 
%         delaysamples = round(Sdelay*Sfs);
%         signalpatch = Ssignal (numel(Ssignal)+1-delaysamples:end);
%         signalpatch = fliplr(signalpatch);
%         Ssignal2 = [Ssignal(delaysamples:end) signalpatch];
%         if numel(Ssignal2)>Ssignal
%             Ssignal2 = Ssignal2(1:numel(Ssignal));
%         end
%         Ssignal = Ssignal2; clear Ssignal2
%     end
%     downsampling = 1;
%     if downsampling 
%         newfs = 500; oldfs = Sfs;
%         [Ssignal,Sfs] = BU_ResampleData(double(Ssignal)',newfs,oldfs,params.recsystem);
%         clear newfs oldfs
%     end
% 
%     ecog.smartlab.tdt.co2 = single(Ssignal); % collect
% end
% save('/home/david/Documents/RRET/patients/NS217/vccp/mat/VCCP_NS217_B30_P1_ecog_epoch.mat', 'ecog');