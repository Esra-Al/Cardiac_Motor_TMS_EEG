clc; clear;
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/';
savepath = '/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/events/'; mkdir(savepath); cd(savepath)
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
conds={'motor'}; bl=1;
%onset_ind=NaN(subnum,30); onset_time=NaN(subnum,30);
load /data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/onsets
%% 
for s = 1:subnum
    subid=subj_names{s};
    EEG=pop_loadset([eegpath subj_names{s} '_' conds{bl} '.set']);
    [EEGout,ind]=pop_selectevent(EEG,'type',{'104'});
    visual_event=cell2mat({EEG.event(ind).latency})/EEG.srate; %in secs
    motor_event=visual_event+onset_time(s,:)/1000; %in secs
    save([savepath '/' 'motor_event' subid], 'motor_event')
end

%%
%     for i=1:length(ind)
%         n_events=length(EEG.event);
%         EEG.event(n_events+1).type='motor_onset';
%         EEG.event(n_events+1).latency=motor_event(i)*EEG.srate;
%     end
%     EEG=eeg_checkset(EEG,'eventconsistency');
%     
%     event_phase=EEG.event;
%     save([save2folder '/' 'EEG_event_phase' subid], 'event_phase')