clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
eeglab
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGrest/';
eegfolder = [eegpath 'merge'];  cd(eegfolder);
save2folder= [eegfolder '/HEP_preprocess']; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
cfg2=[]; cfg2.baseline=[-0.1 0];

%%
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_rest_ica_clean.set']);
    
    EEG=pop_chanedit(EEG, 'append',63,'changefield',{64 'labels' 'M1'},'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG=pop_chanedit(EEG, 'setref',{'1:64' 'M1'});
    EEG = eeg_checkset( EEG );
    
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'M1'},'type',{''},'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',{118.5549},'urchan',{64},'ref',{'M1'},'datachan',{0}));
    EEG = eeg_checkset( EEG );
    EEG = pop_reref( EEG, [31 64] ); %m2 and m1
    EEG.icachansind=[];
   
    %add ECG electrode back
    ECG = pop_loadset([subfolder subid '_tmsall_rest.set']);
    ECG= pop_select(ECG, 'channel', {'ECG'});
    EEG.data(end+1,:,:) = ECG.data;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'ECG';
    EEG = eeg_checkset(EEG);
    
    % Downsample data (5000 Hz to 500 Hz)
    EEG = pop_resample( EEG, 500);
    [EEG, keep_ind] = pop_epoch(EEG, {'ECG'} , [-0.1 0.6]);
    
    data = eeglab2fieldtrip(EEG, 'timelockanalysis','chan_loc');
    data.dimord='chan_time';
    data=ft_timelockbaseline(cfg2, data);
    data.fsample=EEG.srate;
    enum=size(EEG.data,3);
    grandtrials{s,1}=data;
    epochnum(s,1)=enum;
    
    clear EEG data enum subdata2 delind delst keep_ind hep1 hep3 inds
end
save('grandtrials_mastoid_downsample','grandtrials','subnum','epochnum')
%%
clc; clear;
save2folder='/data/p_02186/TMS_ECG2/analyses/EEGrest/merge/HEP_preprocess'; cd(save2folder)
load('grandtrials_mastoid_downsample')
hep=grandtrials(:,1)';
cfg = []; cfg.channel   = 'all'; cfg.latency   = 'all'; cfg.parameter = 'avg';  %cfg.keepindividual='yes';
GA_hep = ft_timelockgrandaverage(cfg,hep{:});
%% figure
load('/data/pt_02333/matlab_code/clusterch_late.mat')
cfg = [];
cfg.channel=chans;
cfg.xlim =[-0.1 0.6]; 
%cfg.ylim =[-9 3];
cfg.ylim =[-2 2];
cfg.linewidth=2;

ft_singleplotER(cfg, GA_hep);
xlabel('Time (s)', 'FontSize', 15);
ylabel('HEP Amplitude (uV)', 'FontSize', 15);
title('');
set(gcf, 'Color', 'w')
ax1=gca;
%set(ax1,'XTick',-1:0.2:2,'FontSize',15)
%set(ax1,'YTick',-1.2:1.2:1.2,'FontSize',15)
%%
 load('/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/HEP_fdi/layout.mat')
cfg = [];
cfg.layout=layout;
cfg.xlim =[-0.05 0.4]; %cfg.ylim =[-0.5 0.5];
cfg.linewidth=2;

ft_multiplotER(cfg, GA_hep);
%xlabel('Time (s)', 'FontSize', 25); ylabel('Potential (�V)', 'FontSize', 25);
title('');
ax1=gca;