clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
addpath('/data/hu_esraal/Documents/TMSEEG-4.0')
addpath('/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/')
eeglab
%%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/old_analyses_tmseeg';  cd(eegfolder);
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/HEP_preTMS_ica2/tmsartefact'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
cfg2=[]; cfg2.baseline=[-0.1 0];
load('events')
%%
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_default_1_2_3_4_5_6_7_8.set']);
    EEG = tmseeg_addTMSTimeBack_interp(EEG, EEG.epoch_length);
    
    %rereference to mastoids
    EEG=pop_chanedit(EEG, 'append',63,'changefield',{64 'labels' 'M1'},'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG=pop_chanedit(EEG, 'setref',{'1:64' 'M1'});
    EEG = eeg_checkset( EEG );
    
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'M1'},'type',{''},'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',{118.5549},'urchan',{64},'ref',{'M1'},'datachan',{0}));
    EEG = eeg_checkset( EEG );
    EEG = pop_reref( EEG, [31 64] ); %m2 and m1
    EEG.icachansind=[];
   
    %load ECG electrode 
    ECG = pop_loadset([subfolder subid '_tmsall_default.set']);
    ECG = pop_epoch(ECG, {'TMS'} , [-1.4 1]);
    ECG = pop_rmbase(ECG, [-110 -10]);
    ECG= pop_select(ECG, 'channel', {'ECG'});
    
    %add ECG electrode back in the EEG structure
    EEG.data(end+1,:,:) = ECG.data;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'ECG';
    EEG = eeg_checkset(EEG);
    
    EEG.event = events{s};
    % Downsample data (5000 Hz to 500 Hz)
    EEG = pop_resample( EEG, 500);
    [EEG, keep_ind] = pop_epoch(EEG, {'fake_ecg'} , [-0.1 0.4]);
    EEG.data = mean(EEG.data,3);
    EEG.trials = 1;
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsartefact.set'], 'filepath', save2folder, 'version', '7.3');

    
    data = eeglab2fieldtrip(EEG, 'timelockanalysis','chan_loc');
    data.dimord='chan_time';
    data=ft_timelockbaseline(cfg2, data);
    data.fsample=EEG.srate;
    enum=size(EEG.data,3);
    grandtrials{s,1}=data;
    epochnum(s,1)=enum;
    
    clear EEG data enum
end
save('grandtrials_tms','grandtrials','subnum','epochnum');
%%
clc; clear;
load('/data/pt_02333/matlab_code/clusterch_late.mat')
cfg = []; cfg.channel   = 'all'; cfg.latency   = 'all'; cfg.parameter = 'avg';  %cfg.keepindividual='yes';
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/HEP_preTMS_ica2/tmsartefact'; cd(save2folder);

load('grandtrials_tms')
hep_tms=grandtrials(:,1)';
GA_hep_tms= ft_timelockgrandaverage(cfg,hep_tms{:});
%%
cfg2 = [];
cfg2.channel=chans;
cfg2.xlim =[-0.1 0.6]; 
cfg2.ylim =[-3 3];
cfg2.linewidth=2;
cfg2.linecolor=[rgb('blue'); rgb('red'); rgb('darkgreen'); rgb('purple'); rgb('pink'); rgb('lightblue')];
ft_singleplotER(cfg2, GA_hep_tms);

set(gcf, 'Color','w')
xlabel('Time (s)'); 
ylabel('HEP Potential (uV)');
%legend('tms','tms_05', 'Location', 'bestoutside', 'box', 'off' )