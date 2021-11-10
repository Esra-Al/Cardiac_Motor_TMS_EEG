clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/hu_esraal/Documents/TMSEEG-4.0')
eeglab
%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGrest/merge/'; cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
%%
s=36
subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
ECG = pop_loadset([subfolder subid '_tmsall_rest.set']);
ECG= pop_select(ECG, 'channel', {'ECG'});

EEG = pop_loadset([subfolder subid '_tmsall_rest_ica.set']);

EEGorig=EEG;

EEG.data(end+1,:,:) = ECG.data;
EEG.nbchan = size(EEG.data,1);
EEG.chanlocs(end+1).labels = 'ECG';
EEG = eeg_checkset(EEG);
EEG = pop_epoch(EEG, {'ECG'} , [-0.05 0.4]);
EEG=pop_rmbase(EEG, [-0.05 0]);


% use SASICA extension of eeglab
close all
cfg.opts.noplot = 0;
cfg.opts.nocompute = 0;
cfg.autocorr.enable = true;
cfg.autocorr.dropautocorr = 'auto'; % default is 2SD: 'auto'
cfg.autocorr.autocorrint = 20;% will compute autocorrelation with this many milliseconds lag
cfg.focalcomp.enable = true;
cfg.focalcomp.focalICAout = 'auto 1.5'; % using default, which is 2SD
cfg.trialfoc.enable = false;
cfg.trialfoc.focaltrialout = 'auto'; % using default which is 2SD
cfg.resvar.enable = false;
cfg.resvar.thresh = 40 ;% %residual variance allowed, default is 15%
cfg.SNR.enable = false;

cfg.EOGcorr.enable = true;
cfg.EOGcorr.corthreshV ='auto 4';% threshold correlation with vertical EOG, default would be 4 SD from average correlation: 'auto 4'
cfg.EOGcorr.Veogchannames = [1];% vertical channel(s), VEOG Fp1     EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_default_1_2_3_4_5_6_7_8_9_10.set'], 'filepath', subfolder, 'version', '7.3');

cfg.EOGcorr.corthreshH ='auto 4';% threshold correlation with horizontal EOG
cfg.EOGcorr.Heogchannames = [11 12];% horizantal channel(s) f7 f8
cfg.chancorr.enable =true;
cfg.chancorr.channames ='ECG';
cfg.chancorr.corthresh='auto 1.5';

cfg.FASTER.enable = false; % default thresh is 3 SDs from average for each measure
cfg.FASTER.blinkchans = [];
cfg.ADJUST.enable = false; % combines different features for detection (SAD,SVD,TK
cfg.MARA.enable = false;
[EEG, cfg] = eeg_SASICA(EEG,cfg);
% this loop opens every ica component
%EEG=pop_select(EEG, 'time', [-0.1 0.4]);
% icomp = 1;
%%
inds=find(EEG.reject.gcompreject==1);
for icomp =  inds
    eeg_SASICA(EEG, ['pop_prop( EEG, 0, ' num2str(icomp) ', findobj(''tag'',''comp' num2str(icomp) '''), { ''freqrange'', [1 50] })']);
end
% for icomp =  sort(1:length(EEG.reject.gcompreject), 'descend')
%     eeg_SASICA(EEG, ['pop_prop( EEG, 0, ' num2str(icomp) ', findobj(''tag'',''comp' num2str(icomp) '''), { ''freqrange'', [1 50] })']);
% end
%%
EEGorig.comprej=find(EEG.reject.gcompreject == 1);
%  newecg=[34]; % ica component to add
%  EEGorig.comprej=unique(sort([EEGorig.comprej newecg]));

x = inputdlg('Enter ecg components:',...
    'Sample', [1 50]);
EEGorig.ecgrej=str2num(x{:});
disp(EEGorig.comprej)
%
EEGclean = pop_subcomp(EEGorig,EEGorig.comprej, 0);
EEGclean = eeg_checkset(EEGclean);
pop_saveset(EEGclean, 'filename', [subid '_tmsall_rest_ica_clean.set'], 'filepath', subfolder, 'version', '7.3');