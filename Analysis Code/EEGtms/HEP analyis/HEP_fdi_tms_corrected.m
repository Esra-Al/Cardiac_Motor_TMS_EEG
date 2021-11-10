clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
addpath('/data/hu_esraal/Documents/TMSEEG-4.0')
addpath('/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/')
eeglab
%%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/old_analyses_tmseeg';  cd(eegfolder);
save2folder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/HEP_fdi/tms_corrected'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
conds={'mep1','mep2', 'mep3'}; % lowmep=mep1 
condsnum=length(conds);
grandtrials=cell(subnum,condsnum); epochnum=NaN(subnum,condsnum);
cfg2=[]; cfg2.baseline=[-0.1 0];
tmsfolder= '/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/HEP_preTMS_ica2/tmsartefact';
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
    
    % Downsample data (5000 Hz to 500 Hz)
    EEG = pop_resample( EEG, 500);
    
    ecginds=find(strcmp('ECG',{EEG.event.type}));
    ecglat=[EEG.event(ecginds).latency];
    tmslat=[EEG.event(strcmp('TMS',{EEG.event.type})).latency];
    delst=[];
    for st=1:length(tmslat)
        pre_ecg= find(tmslat(st)>ecglat, 1, 'last');
        rind1=ecginds(pre_ecg);
        
        if (tmslat(st)-ecglat(pre_ecg))/EEG.srate >0.4 && (tmslat(st)-ecglat(pre_ecg))/EEG.srate <1.4
            EEG.event(rind1).type= 'pre_ECG';
        else
            delst=[delst st];
        end
    end
    %
    
    [EEG, keep_ind] = pop_epoch(EEG, {'pre_ECG'} , [-0.1 0.4]);
    faketms=pop_loadset([subid '_tmsartefact.set'], tmsfolder);
    EEG.data= bsxfun(@minus, EEG.data, faketms.data);
        
    delind=setdiff(1:416-length(delst), keep_ind);
    
    %
    load(['/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/FDIdata/data' subid '.mat'])
    subdata2(:,3)=1:416;
    subdata2(delst,:)=[]; 
    if delind
        subdata2(delind,:)=[]; 
    end
    subdata2(:,4)=quantileranks(subdata2(:,2),3,0); %3 bins
    
    hep1=find(subdata2(:,4)==1);  hep2=find(subdata2(:,4)==2); hep3=find(subdata2(:,4)==3);
    inds={hep1,hep2, hep3};
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
    clear subdata2 delind delst keep_ind hep1 hep3 inds
end
   
save('grandtrials','grandtrials','subnum','epochnum')
%%
clc; clear;
save2folder='/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/reref/HEP_fdi/tms_corrected'; cd(save2folder)
load('grandtrials')
meplow=grandtrials(:,1)';
mepmed=grandtrials(:,2)';
mephigh=grandtrials(:,3)';
cfg = []; cfg.channel   = 'all'; cfg.latency   = 'all'; cfg.parameter = 'avg';  %cfg.keepindividual='yes';
GA_meplow = ft_timelockgrandaverage(cfg,meplow{:});
GA_mepmed = ft_timelockgrandaverage(cfg,mepmed{:});
GA_mephigh = ft_timelockgrandaverage(cfg,mephigh{:});
load('/data/pt_02333/matlab_code/clusterch_late.mat')
%%
cfg = [];
cfg.elec =meplow{1,1}.elec;
cfg.rotate = 90;
cfg.channel     = {'all',  '-VEOG', '-ECG'}; % exclude EOG channel
layout = ft_prepare_layout(cfg);
layout.pos(:,2) = (layout.pos(:,2)+0.02).*0.92 ;
cfg.layout=layout;
ft_layoutplot(cfg)
save layout.mat layout

cfg.method      = 'triangulation';
cfg.feedback    = 'no'; % don't show a neighbour plot
neighbours= ft_prepare_neighbours(cfg, meplow{1,1}); % define neighbouring channels
cfg.neighbours = neighbours;
save('neighbours.mat', 'neighbours');
%% cluster based permutation
load('layout.mat'); load('neighbours.mat');
cfg = [];
cfg.channel     = chans; % exclude EOG channel
%cfg.channel     = {'all',  '-ECG'}; % exclude EOG channel
cfg.latency = [0.296 0.4];
%cfg.avgovertime = 'yes'; cfg.parameter   = 'avg';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; % within subnum design
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric';
cfg.minnbchan = 2; % at least two adjacent channels to enter consideration
cfg.neighbours = neighbours; % as defined above
cfg.tail = 0; % two-tailed test
cfg.clustertail = 0;
cfg.alpha = 0.025;
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
stat= ft_timelockstatistics(cfg, meplow{:}, mephigh{:});
cfg = [];
cfg.layout = layout;
%cfg.highlightseries={'on','labels'}; %cfg.subplotsize=[1 1];cfg.visible  = 'off';
addpath '/data/hu_esraal/Documents/eeglab2019_0/plugins/Fieldtrip-lite20210209/'
ft_clusterplot_noplot(cfg,stat);
%%
%poschans=stat.label(sum(stat.posclusterslabelmat==1,2)>0);
negchans=stat.label(sum(stat.negclusterslabelmat==1,2)>0);
cfg = [];
cfg.channel=negchans;
cfg.xlim =[-0.1 0.4]; cfg.ylim =[-2 2];
cfg.baseline=[-0.1 0]
cfg.graphcolor='rb';
cfg.linewidth=2;
ft_singleplotER(cfg, GA_meplow, GA_mephigh);
legend('low', 'high')
%% figure
 load('/data/pt_02333/matlab_code/clusterch_late.mat')
cfg = [];
cfg.channel=chans;
cfg.xlim =[-0.05 0.42]; cfg.ylim =[-2 2];
%cfg.baseline=[-0.05 0]
cfg.linecolor=[rgb('pink'); rgb('purple')];
cfg.linewidth=2;
%hFig=figure(1);
%set(hFig, 'Position', [100 100 800 800])
ft_singleplotER(cfg, GA_meplow, GA_mephigh);

legend('low MEP', 'high MEP', 'Location','bestoutside')
xlabel('Time (s)', 'FontSize', 15); 
ylabel('HEP Potential (UV)', 'FontSize', 15); 
title('');

[id, x] = match_str(GA_meplow.label, cfg.channel);
hold on
%
a=[0.304 0.328]-0.001;
xall=GA_meplow.time(GA_meplow.time >=a(1) & GA_meplow.time <=a(2));
y1=mean(GA_meplow.avg(id, GA_meplow.time >=a(1) & GA_meplow.time <=a(2)));
y2=mean(GA_mephigh.avg(id, GA_meplow.time >=a(1) & GA_meplow.time <=a(2)));
yall=[y1, fliplr(y2)];
xall=[xall,sort(xall, 'descend')];
fill(xall,yall,[0.5 0.5 0.5], 'LineStyle','none')
alpha(0.4)
set(gcf,'Color','w');
saveas(gcf, 'hep_meps.svg')
%%
load('layout.mat')
cfg.operation='subtract';
cfg.parameter = 'avg';

GAhigh_low=ft_math(cfg, GA_mephigh, GA_meplow);

cfg=[];
cfg.layout=layout;
cfg.style = 'straight'; cfg.comment= 'no'; cfg.marker = '';
cfg.zlim = [-0.8 0.8]; cmap=colormap("jet"); cfg.colormap=cmap;
%figure('units','normalized','outerposition',[0 0 1 1])
figure
cfg.xlim = [0.304 0.328]; %cfg.colorbar = 'southoutside';
ax1=gca;
set(ax1,'FontSize',15,'Color','w')
ft_topoplotER(cfg,GAhigh_low);
%saveas(gcf,'.svg')
print(gcf, 'hep_topo.png', '-dpng', '-r900');
