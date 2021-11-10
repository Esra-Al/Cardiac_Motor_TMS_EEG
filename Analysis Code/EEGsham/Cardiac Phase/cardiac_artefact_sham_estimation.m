clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
eeglab
%%
ts=[-100 600];     % <--- input in ms
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/';  cd(eegfolder);
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase_fake/'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
cfg2=[]; cfg2.baseline=[-0.1 0];

conds={'fake_sys','fake_dys'};
condsnum=length(conds);
%%
for s=1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_firstica2_interpolate_filter.set']);

    load([save2folder '/' 'EEG_event_phase' subid '.mat']) %load fake trigger event info for each subject
    EEG = pop_resample( EEG, 500);
    EEG.event=event_phase; %replace event structure
    for c=1:condsnum
        hd_epoch_HER(EEG, conds(c),ts, subid, save2folder)
    end
end
% grand average
cd(save2folder)
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
cfg2=[]; cfg2.baseline=[-0.1 0];
for s=1:subnum
    subid = subj_names{s};
    for c=1:condsnum
        [data,enum]=hd_grand([conds{c} '.set'], subid, save2folder,cfg2);
        grandtrials{s,c}=data;
        epochnum(s,c)=enum;
    end
end
save('grandtrials','grandtrials','conds','subnum','epochnum')
%%
clc; clear;
save2folder='/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase_fake/'; 
cd(save2folder)
load('grandtrials')
fake_sys=grandtrials(:,1)';
fake_dys=grandtrials(:,2)';

cfg = []; cfg.channel   = 'all'; cfg.latency   = 'all'; cfg.parameter = 'avg';
GA_fake_sys = ft_timelockgrandaverage(cfg, fake_sys{:});
GA_fake_dys= ft_timelockgrandaverage(cfg,fake_dys{:});
%%
cfg = [];
cfg.channel={'C4'}; cfg.xlim =[-0.1 0.6]; cfg.ylim=[-0.5 1];
cfg.graphcolor='rb';
cfg.linewidth=2;
figure; ft_singleplotER(cfg, GA_fake_sys,GA_fake_dys);
xlabel('Time (s)', 'FontSize', 25); ylabel('Potential (�V)', 'FontSize',25);
title('');
ax1=gca;
%set(ax1,'XTickLabel', '', 'YTickLabel', '','FontSize',15)
legend('systole', 'diastole')
set(gcf,'Color','w');
%export_fig('poster_fakesys_dys_C4')
%saveas(gcf, 'allsysdysgrand_C4.svg')