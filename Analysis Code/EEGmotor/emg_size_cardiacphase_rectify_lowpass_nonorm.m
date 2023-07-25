clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
eeglab
%%
ts=[-1000 4000];
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/';
motorpath='/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/events/';
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/merge/motor_onset/emg_cardiacphase/'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
conds={'sys','dys','all'};
condsnum=length(conds);
load("/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/sys_dys_motor_onset.mat")
sys_dys=x;
grandtrials=cell(subnum,condsnum); epochnum=NaN(subnum,condsnum);
cfg2=[]; cfg2.baseline=[-0.11 -0.01];
% add trial numbers for each subject
for s=1:subnum
    sys_dys(sys_dys(:,2)==s,1)=1:length(sys_dys(sys_dys(:,2)==s,1));
end
%%
for s=1:subnum
    subid=subj_names{s};
    EEG=pop_loadset([eegpath subj_names{s} '_motor.set']);
    EEG=pop_select(EEG, 'channel', {'FDI'});
    load([motorpath '/' 'motor_event' subid])
    
    for i=1:length(motor_event)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='motor_onset';
        EEG.event(n_events+1).latency=motor_event(i)*EEG.srate;
    end
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    sig_sys=sys_dys(sys_dys(:,2)==s & sys_dys(:,5)==1,1); % this code is different than tms-eeg !!
    sig_dys=sys_dys(sys_dys(:,2)==s & sys_dys(:,6)==1,1);
    inds={sig_sys, sig_dys, 1:length(motor_event)};
    
    [b,a]=butter(2,[10, 500]/(EEG.srate/2)); % 10 500
    EEG.data=filtfilt(b,a,double(EEG.data)')';
    [c,d]=butter(2,[45, 55]/(EEG.srate/2),'stop'); % notch
    EEG.data=filtfilt(c,d,double(EEG.data)')';
    EEG.data=abs(EEG.data);
    
    [e,f]  = butter(4, 8/(EEG.srate/2),'low'); 
    EEG.data=filtfilt(e,f,double(EEG.data)')';
    EEG=pop_resample(EEG,500);
    EEG=pop_epoch(EEG,{'motor_onset'},ts/1000,'epochinfo','yes'); %try with baseline after this!
    EEG=pop_rmbase(EEG, [-110 -10]);
    for c=1:condsnum
        EEG_cond=pop_select(EEG, 'trial', inds{c});
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
save('grandtrials_linenv_4th_8hz_base','grandtrials','conds','subnum','epochnum')
%%
clc; clear;
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/merge/motor_onset/emg_cardiacphase/'; cd(save2folder)
load('grandtrials_linenv_4th_8hz_base')
subnum=length(grandtrials);
sig_sys=grandtrials(:,1)';
sig_dys=grandtrials(:,2)';
sig_all=grandtrials(:,3)';
cfg = []; cfg.channel   = 'all'; cfg.latency   = 'all'; cfg.parameter = 'avg';
GA_sys = ft_timelockgrandaverage(cfg, sig_sys{:});
GA_dys = ft_timelockgrandaverage(cfg,sig_dys{:});
GA_all = ft_timelockgrandaverage(cfg,sig_all{:});
load('/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase/layout.mat')
load('/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase/neighbours.mat')
%%
cfg = [];
cfg.channel={'FDI'};
cfg.xlim =[-1 3]; %cfg.ylim=[-4 7];
cfg.graphcolor='rbg';
cfg.linewidth=2;
ft_singleplotER(cfg, GA_sys_diff,GA_dys_diff); legend('systole', 'diastole')
%%
systole=NaN(1,subnum); distole=NaN(1,subnum);
inds=find(sig_sys{1}.time >=0.34 & sig_sys{1}.time <=0.454);

for s=1:subnum
    systole(s)=mean(sig_sys{s}.avg(inds));
    diastole(s)=mean(sig_dys{s}.avg(inds));
end
%[h,p]=ttest(systole, diastole)
cohensd = mean(systole-diastole) ./ std(systole-diastole)
%%
cfg = [];
cfg.channel={'FDI'};
cfg.xlim =[-1 4]; %cfg.ylim=[-4 7];
cfg.graphcolor='rbg';
cfg.linewidth=2;
% ft_singleplotER(cfg, GA_sys,GA_dys, GA_all); legend('systole', 'diastole', 'all')
ft_singleplotER(cfg, GA_sys,GA_dys); %legend('systole', 'diastole')
%xlabel('Time (s)', 'FontSize', 25); %ylabel('Potential (�V)', 'FontSize', 25);
title(' ');
ax1=gca;
% set(ax1,'XTick', -0.1:0.1:0.6,'FontSize',12,'Color','w')
hold on
a=[0.34 0.454]-0.001;
xall=GA_sys.time(GA_sys.time >=a(1) & GA_sys.time <=a(2));
y2=GA_sys.avg(GA_sys.time >=a(1) & GA_sys.time <=a(2));subnum=length(grandtrials);
y1=GA_dys.avg(GA_sys.time >=a(1) & GA_sys.time <=a(2));
yall=[y1, fliplr(y2)];
xall=[xall,sort(xall, 'descend')];
fill(xall,yall,[0.5 0.5 0.5], 'LineStyle','none')
alpha(0.4)
set(gcf,'Color','w');
saveas(gcf,'motoronset_emg_sys_dys_revision.svg')
%%
cfg = [];
%cfg.channel     = {'all', '-ECG'}; % exclude ECG channel
cfg.channel={'FDI'};
cfg.latency = [0 1]; %
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; % within subnum design
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic ='maxsum';
cfg.clusterthreshold = 'parametric';
cfg.neighbours = neighbours; % as defined above
cfg.tail = 0; % two-tailed test
cfg.clustertail = 0;
cfg.numrandomization = 1000;
% prepare design matrix for comparison of two conditions

design = zeros(2,2*subnum);
for i = 1:subnum
    design(1,i) = i;
end
for i = 1:subnum
    design(1,subnum+i) = i;
end
design(2,1:subnum)        = 1;
design(2,subnum+1:2*subnum) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
% Run statistics

stat= ft_timelockstatistics(cfg, sig_sys{:}, sig_dys{:});
cfg = [];
cfg.alpha=0.1;
cfg.layout = layout;
ft_clusterplot_noplot(cfg,stat);
