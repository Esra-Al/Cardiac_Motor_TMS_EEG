clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
eeglab
%%
ts=[-100 600];
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/old_analyses_tmseeg';  cd(eegfolder);
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/cardiacphase_fake/'; mkdir(save2folder); cd(save2folder);
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
    
%     EEG=pop_chanedit(EEG, 'append',63,'changefield',{64 'labels' 'M1'},'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
%     EEG=pop_chanedit(EEG, 'setref',{'1:64' 'M1'});
%     EEG = eeg_checkset( EEG );
%     
%     EEG = pop_reref( EEG, [],'refloc',struct('labels',{'M1'},'type',{''},'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',{118.5549},'urchan',{64},'ref',{'M1'},'datachan',{0}));
%     EEG = eeg_checkset( EEG );
%     EEG = pop_reref( EEG, [31 64] ); %m2 and m1
%     
    EEG = pop_resample(EEG, 500);
    EEG.event=event_phase; %replace event structure
     
    for c=1:condsnum
        hd_epoch_HER(EEG, conds(c),ts, subid, save2folder)
    end
end
% grand average
cd(save2folder)
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
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
save2folder='/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/cardiacphase_fake/'; 
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