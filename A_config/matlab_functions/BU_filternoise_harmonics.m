function [eeg] = BU_filternoise_harmonics(eeg,fs,filternoise_harmonics,selchan,plotting)
%%IN:
%eeg ; chan*samples ; continous data
%fs = 512
%filternoise_harmonics=[60 120 180];%[60];
%selchan = 1; %%chan 2 plot

% addpath('/Users/jherreroru/Documents/gitmatlab/fct/');
eeg = fillmissing(eeg, 'linear', 2);   %sam added to handle NaNs in eeg

if isempty(selchan)
    selchan = 1:size(eeg,1);
end

%% -------------------------------------------------- 
%% visualize raw signal (select chan) 
if plotting
    close all;
    curfig = figure;
    set(curfig,'position',[100   50   1500   1200],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
    x = eeg(selchan,:);
   
    %XX = PSD(x,fs,257,round(length(eeg)/257),20);
    [pxx,hzPxx] = pwelch(x,4*fs,2*fs,4*fs,fs);
    %BU_plot_1chan_FFT (eeg(2,:),ecog.eeg.fs,length(eeg(1,:))/ecog.eeg.fs  );
    subplot(221);plot(x); axis tight; hold on;
    subplot(223);plot(x); axis tight;xlim([fs fs*10]);hold on; %10 sec
    subplot(254);semilogy(hzPxx,pxx);xlim([0 75]);hold on;subplot(255);plot(hzPxx,pxx);xlim([45 70]); hold on;
    disp(class(x)); disp(size(x)); whos x
    %plotpsd(x,fs,'logx')
    
end


%% -------------------------------------------------- 
%% filter noise 

eegF=eeg;
if ~isempty(filternoise_harmonics)
    %filternoise_harmonics=[60:60:180];
    for i=1:length(filternoise_harmonics)
        wo = filternoise_harmonics(i)/(fs/2);
        bw = wo/50;%Q-factor=100
        [b,a] = iirnotch(wo,bw,20);
        for chan=1:size(eeg,1)
            X = double(eegF(chan,:));
            eegF(chan,:)=filtfilt(b,a,X);
        end
    end
end
eeg = eegF;clear eegF X b a wo;
%% -------------------------------------------------- 
%% visualize notch-free signal (select chan) 
plotting=1;
if plotting
    %close all hidden
    %curfig = figure;
    %set(curfig,'position',[100   50   1200   400],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
    x = eeg(selchan,:);
    %XX = PSD(x,fs,257,round(length(eeg)/257),20);
    [pxx,hzPxx] = pwelch(x,4*fs,2*fs,4*fs,fs);
    subplot(221);plot(x,'r');hold on; axis tight
    subplot(223);plot(x,'r');hold on; axis tight;xlim([fs fs*10]);%10 sec
    subplot(254);semilogy(hzPxx,pxx);xlim([0 75]);hold on;subplot(255);plot(hzPxx,pxx);xlim([45 70]); hold on;
end


rmpath('/Users/jherreroru/Documents/gitmatlab/fct/');