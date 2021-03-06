clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
addpath('/data/p_02186/TMS_ECG2/analyses/EEGmotor/merge/cardiacphase_fakecorr/time_frequency')
eeglab
%%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/merge/';
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/merge/motor_onset/cardiacphase_fakecorr/time_frequency'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
cfg2=[]; cfg2.baseline=[-0.1 0];

conds={'sys','dys'};
condsnum=length(conds);
load("/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/sys_dys_motor_onset.mat")
sys_dys=x;
grandtrials=cell(subnum,condsnum); epochnum=NaN(subnum,condsnum);
cfg2=[]; cfg2.baseline=[-0.1 0];
%fakefolder='/data/p_02186/TMS_ECG2/analyses/EEGrest/merge/cardiacphase_fake_tms_events3/';
% add trial numbers for each subject
for s=1:subnum
    sys_dys(sys_dys(:,2)==s,1)=1:length(sys_dys(sys_dys(:,2)==s,1));
end
%%
for s=1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_motor_ica_clean.set']);
   
    EEG=pop_chanedit(EEG, 'append',63,'changefield',{64 'labels' 'M1'},'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG=pop_chanedit(EEG, 'setref',{'1:64' 'M1'});
    EEG = eeg_checkset( EEG );
    
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'M1'},'type',{''},'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',{118.5549},'urchan',{64},'ref',{'M1'},'datachan',{0}));
    EEG = eeg_checkset( EEG );
    EEG = pop_reref( EEG, [31 64] ); %m2 and m1
    
    sig_sys=sys_dys(sys_dys(:,2)==s & sys_dys(:,5)==1,1); % this code is different than tms-eeg !!
    sig_dys=sys_dys(sys_dys(:,2)==s & sys_dys(:,6)==1,1);
    %numtrial=min([length(sig_sys) length(sig_dys)]); sig_sys=sig_sys(1:numtrial); sig_dys=sig_dys(1:numtrial);
    inds={sig_sys, sig_dys};
    %: Downsample data (5000 Hz to 500 Hz)
    EEG = pop_resample( EEG, 500);
    [EEG, keep_ind] = pop_epoch(EEG, {'104'} , [-1 3]);
    for c=1:condsnum 
        EEG_cond=pop_select(EEG, 'trial', inds{c});
        [TFA]=eegl_calcTFA_3(EEG_cond,5:40, [4 10], 1, 'electrode');
        save(['TFA_' subid conds{c} '.mat'], '-struct', 'TFA');  
        
        data = eeglab2fieldtrip(EEG_cond, 'timelockanalysis','chan_loc');
        data.dimord='chan_time';
        data=ft_timelockbaseline(cfg2, data); 
        data.fsample=EEG_cond.srate;
        enum=size(EEG_cond.data,3);
        grandtrials{s,c}=data;
        epochnum(s,c)=enum;  
        clear EEG_cond data enum 
    end   
    
end
save('grandtrials','grandtrials','conds','subnum','epochnum', 'subj_names')
%%
clc; clear;
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/merge/motor_onset/cardiacphase_fakecorr/time_frequency'; cd(save2folder);
load('grandtrials.mat')

label=grandtrials{1,1}.label; 
TFsys=field_tf('sys', subj_names, label);
TFdys=field_tf('dys', subj_names, label);
cfg=[];  cfg.baseline=[-0.110 -0.01]; cfg.baselinetype='absolute';
TFsys=ft_freqbaseline(cfg, TFsys); TFdys=ft_freqbaseline(cfg, TFdys); 

cfg=[]; cfg.operation='subtract'; cfg.parameter='powspctrm';
TFsys_dys=ft_math(cfg, TFsys, TFdys);

load('/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase/layout.mat')
load('/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase/neighbours.mat')
%%
grandsys=TFsys;
grandsys.powspctrm=squeeze(mean(TFsys.powspctrm,1));
grandsys.dimord='chan_freq_time' ;

granddys=TFdys;
granddys.powspctrm=squeeze(mean(TFdys.powspctrm,1));
granddys.dimord='chan_freq_time' ;

grandsys_dys=TFsys_dys;
grandsys_dys.powspctrm=squeeze(mean(TFsys_dys.powspctrm,1));
grandsys_dys.dimord='chan_freq_time' ;
%%
cfg = [];
cfg.zlim = [-1 1];	    
cfg.xlim= [-0.1 1];
cfg.channel={'C4';'CP6';'CP4';'C6'};
cmap=colormap("jet"); cfg.colormap=cmap;
% ft_singleplotTFR(cfg, grandsys);
% ft_singleplotTFR(cfg, granddys);
ft_singleplotTFR(cfg, grandsys_dys);
saveas(gcf,'alpha_emg_sys_dys.svg')
%% Cluster stats
cfg = []; cfg.latency = [0 1];
cfg.channel={'C4';'CP6';'CP4';'C6'};
cfg.frequency = [8 30]; %cfg.avgoverfreq='yes';
cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_depsamplesT'; % within subnum design
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum'; 
cfg.neighbours = neighbours; % as defined above
cfg.tail = 0; % two-tailed test
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 1000;

design = zeros(2,2*subnum); % prepare design matrix for comparison of two conditions
for i = 1:subnum
    design(1,i) = i;
end
for i = 1:subnum
    design(1,subnum+i) = i;
end
design(2,1:subnum)        = 1;
design(2,subnum+1:2*subnum) = 2;
cfg.design = design; cfg.uvar  = 1; cfg.ivar  = 2;

[stat] = ft_freqstatistics(cfg, TFsys, TFdys);
%%
cfg = []; 
cfg.highlightseries={'on','labels'};
cfg.parameter = 'stat';
cfg.layout = layout;
%cfg.subplotsize=[1 1];
ft_clusterplot_noplot(cfg, stat);
%%
negchans=stat.label(sum(sum(stat.negclusterslabelmat==1,2),3)>0);
negtime=stat.time(sum(sum(stat.negclusterslabelmat==1,1),2)>0)
negfreq=stat.freq(sum(sum(stat.negclusterslabelmat==1,1),3)>0)
%%
cfg=[];
cfg.layout=layout;
cfg.style = 'straight'; cfg.comment= 'no'; cfg.marker = '';
cmap=colormap("jet"); cfg.colormap=cmap;

% figure('units','normalized','outerposition',[0 0 1 1])
cfg.xlim = [0 0.726]; cfg.colorbar = 'no';
% cfg.ylim = [8 13]; %[8 13]; %14 25
cfg.ylim = [14 25]; 
cfg.zlim = [-0.5 0.5]; 
ft_topoplotTFR(cfg,grandsys_dys);
set(gcf,'Color','w');
%saveas(gcf,'toposys_dys2.svg')


% cb = colorbar('southoutside');
% cb.XTick=-1:0.5:1;
% set(cb,'position',[.28 .01 .5 .08], 'XTickLabel',{-1 0 1} )
% print(gcf, 'topo_alpha_sys_dys_colorbar.png', '-dpng', '-r900');
print(gcf, 'topo_beta_sys_dys_colorbar.png', '-dpng', '-r900');