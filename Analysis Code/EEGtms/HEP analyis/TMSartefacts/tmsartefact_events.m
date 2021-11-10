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
events=cell(subnum,1); 
removeFrom  = -2;
removeTo = 8;
%%
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_default.set']);
    
    ecginds=find(strcmp('ECG',{EEG.event.type}));
    ecglat=[EEG.event(ecginds).latency];
    tmslat=[EEG.event(strcmp('TMS',{EEG.event.type})).latency];
    tmsdistances=NaN(length(tmslat),1);  delst=[];
    for st=1:length(tmslat)
        pre_ecg= find(tmslat(st)>ecglat, 1, 'last');
        rind1=ecginds(pre_ecg);
        if (tmslat(st)-ecglat(pre_ecg))/EEG.srate >0.4 && (tmslat(st)-ecglat(pre_ecg))/EEG.srate <1.4
            tmsdistances(st)=tmslat(st)-ecglat(pre_ecg);
        else
            delst=[delst st];
        end
    end
    tmsdistances(delst)=[];
    tmslat_keep=tmslat; tmslat_keep(delst)=[];
    
    for sh=1:10 %shuffle 10 times
        distance_sh=shuffle(tmsdistances);
        for st=1:length(distance_sh)
            n_events=length(EEG.event);
            EEG.event(n_events+1).type='fake_ecg';
            EEG.event(n_events+1).latency=tmslat(st)-distance_sh(st);
            EEG.event(n_events+1).urevent=n_events+1;
            
        end
    end
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    EEG= pop_select(EEG, 'nochannel', {'ECG','FDI'});
    EEG = pop_epoch(EEG, {'TMS'} , [-1.4 1]);
    EEG = pop_rmbase(EEG, [-110 -10]); 
    
%     temp_removeFrom = find(EEG.times<=removeFrom);
%     temp_removeTo = find(EEG.times<=removeTo);
%     period2remove = temp_removeFrom(end)+1:temp_removeTo(end)-1;
%     EEG = eeg_checkset(EEG);
%     
%     EEG.data(:,period2remove,:) = [];
%     EEG.times(period2remove)    = [];
%     EEG.TMS_period2remove_1     = period2remove;
%     EEG.custom.times            = EEG.times;
%     EEG.pnts                    = numel(EEG.times);
%     EEG                         = eeg_checkset( EEG ); % some adjustments would be made
%     nevents = length(EEG.event);
%     for index = 1 : nevents
%         EEG.event(index).latency=EEG.event(index).latency-(index-1)*length(period2remove);
%     end
    
    events{s}=EEG.event;
    
    clear EEG data enum
end
save('events','events','subnum')
%%
