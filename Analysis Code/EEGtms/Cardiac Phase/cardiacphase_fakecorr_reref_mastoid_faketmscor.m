clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
eeglab
%%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/old_analyses_tmseeg/';  cd(eegfolder);
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/cardiacphase_fakecorr/'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
conds={'systole','diastole'};
condsnum=length(conds);
load("/data/p_02186/TMS_ECG2/analyses/EKGdata/det_diff2peak/ecgdata.mat")
sys_dys=x;
grandtrials=cell(subnum,condsnum); epochnum=NaN(subnum,condsnum);
cfg2=[]; cfg2.baseline=[-0.110 -0.010];
fakefolder='/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/cardiacphase_fake/';
%%
for s=1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_firstica2_interpolate_filter.set']);
    
    EEG=pop_chanedit(EEG, 'append',63,'changefield',{64 'labels' 'M1'},'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG=pop_chanedit(EEG, 'setref',{'1:64' 'M1'});
    EEG = eeg_checkset( EEG );
    
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'M1'},'type',{''},'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',{118.5549},'urchan',{64},'ref',{'M1'},'datachan',{0}));
    EEG = eeg_checkset( EEG );
    EEG = pop_reref( EEG, [31 64] ); %m2 and m1
    
    
    sig_sys=sys_dys(sys_dys(:,2)==s & sys_dys(:,6)==1,1)-(416*(s-1));
    sig_dys=sys_dys(sys_dys(:,2)==s & sys_dys(:,7)==1,1)-(416*(s-1));
    inds={sig_sys, sig_dys};
    %: Downsample data (5000 Hz to 500 Hz)
    EEG = pop_resample( EEG, 500);
    
    for c=1:condsnum 
        EEG_cond=pop_select(EEG, 'trial', inds{c},'time', [-0.1 0.6]);
        ph=conds{c}; ph=ph(end-2:end);
        fakenull=pop_loadset(['VP' subid 'fake_' ph '.set'], fakefolder);

%         % rereference of artefact estimation
        fakenull=pop_chanedit(fakenull, 'append',63,'changefield',{64 'labels' 'M1'},'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
        fakenull=pop_chanedit(fakenull, 'setref',{'1:64' 'M1'});
        fakenull = eeg_checkset( fakenull );
        
        fakenull = pop_reref( fakenull, [],'refloc',struct('labels',{'M1'},'type',{''},'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',{118.5549},'urchan',{64},'ref',{'M1'},'datachan',{0}));
        fakenull = eeg_checkset( fakenull );
        fakenull = pop_reref(fakenull, [31 64] ); %m2 and m1
        
        EEG_cond.data= bsxfun(@minus, EEG_cond.data, mean(fakenull.data,3));
        
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
save('grandtrials_mastoid_base110_faketmscor','grandtrials','conds','subnum','epochnum')
%%
clc; clear;
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/cardiacphase_fakecorr/'; cd(save2folder);
load('grandtrials_mastoid_base110_faketmscor')
subnum=length(grandtrials);
sig_sys=grandtrials(:,1)';
sig_dys=grandtrials(:,2)';
cfg = []; cfg.channel   = 'all'; cfg.latency   = 'all'; cfg.parameter = 'avg';
GA_sys = ft_timelockgrandaverage(cfg, sig_sys{:});
GA_dys = ft_timelockgrandaverage(cfg,sig_dys{:});
load('/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase/layout.mat')
load('/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase/neighbours.mat')
%%
cfg = [];
cfg.channel={'C4';'CP6';'CP4';'C6'};
cfg.xlim =[-0.1 0.3]; cfg.ylim=[-7 7]; 
cfg.graphcolor='rb';
cfg.linewidth=2;
ft_singleplotER(cfg, GA_sys,GA_dys);
%xlabel('Time (s)'); ylabel('Potential (uV)'); 
title(' ');
ax1=gca;
set(ax1,'XTick', -0.1:0.1:0.6,'FontSize',12,'Color','w')
%hline(0,'k'); vline(0,'k')
hold on
%legend('systole', 'diastole')
[id, x] = match_str(GA_sys.label, cfg.channel);

a=[0.022 0.062]-0.001;
xall=GA_sys.time(GA_sys.time >=a(1) & GA_sys.time <=a(2));
y2=mean(GA_sys.avg(id, GA_sys.time >=a(1) & GA_sys.time <=a(2)));
y1=mean(GA_dys.avg(id, GA_sys.time >=a(1) & GA_sys.time <=a(2)));
yall=[y1, fliplr(y2)];
xall=[xall,sort(xall, 'descend')];
fill(xall,yall,[0.5 0.5 0.5], 'LineStyle','none')
alpha(0.4)
set(gcf,'Color','w');
saveas(gcf,'TEPssys_dys2.svg')
%% cluster based permutation
cfg = [];
cfg.channel={'C4';'CP6';'CP4';'C6'};
cfg.latency = [0.015 0.06]; %0.0232 0.060
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
cfg.layout = layout;
ft_clusterplot_noplot(cfg,stat);
close all;
%%

cfg.operation='subtract';
cfg.parameter = 'avg';

GAsys_dys=ft_math(cfg,GA_sys, GA_dys);

cfg=[];
cfg.layout=layout;
cfg.style = 'straight'; cfg.comment= 'no'; cfg.marker = '';
cmap=colormap("jet"); cfg.colormap=cmap;
cfg.zlim = [-1 1]; 
%figure('units','normalized','outerposition',[0 0 1 1])
cfg.xlim = [0.022 0.060]; %cfg.colorbar = 'southoutside';

ft_topoplotER(cfg,GAsys_dys);
set(gcf,'Color','w');
%saveas(gcf,'toposys_dys2.svg')
print(gcf, 'topo_tms_sys_dys_colorbar2.png', '-dpng', '-r900');

cb = colorbar('southoutside');
cb.XTick=-1:0.5:1;
set(cb,'position',[.28 .01 .5 .08], 'XTickLabel',{-1 0 1} )

%%
