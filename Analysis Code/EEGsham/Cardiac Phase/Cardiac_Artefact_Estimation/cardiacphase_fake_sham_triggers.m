clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
eeglab
%%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGsham/merge';
save2folder = '/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase_fake'; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
%% Randomization of the arbitrary triggers
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_sham.set']);
    EEG= pop_select(EEG, 'nochannel', {'ECG','FDI'});
    
    ecginds=find(strcmp('ECG',{EEG.event.type}));
    ecglat=[EEG.event(ecginds).latency];
    
    [EEGout,indb]=pop_selectevent(EEG,'type','boundary');
    [EEGout,indr]=pop_selectevent(EEG,'type','ECG');
    fakelat=NaN(8*(length(indr)-1),1);
    
    for i=1:length(indr)-1
        if sum(indb==indr(i)) || sum(indb==indr(i+1))
            continue;
        else
            for k=sort(0:7, 'descend')
                fakelat(8*i-k)= EEG.event(indr(i)).latency + (EEG.event(indr(i+1)).latency-EEG.event(indr(i)).latency)*rand(1,1);
            end
        end
    end
    fakelat(isnan(fakelat)) = [];
    fake_event=fakelat/EEG.srate; %in secs
    save([save2folder '/' 'fake_event' subid], 'fake_event')
    
    ecg_event=ecglat/EEG.srate; %in secs
    save([save2folder '/' 'ecg_event' subid], 'ecg_event')
end

%% filter ECG
addpath('/data/p_02186/TMS_ECG2/analyses/TESA1.1.1')
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_sham.set']);
    
    EEG=pop_select(EEG, 'channel', {'ECG'});
    EEG = pop_tesa_removedata(EEG, [-2 10],[], {'TMS'} );
    EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );
    
    [c,d]=butter(2,0.5/(EEG.srate/2),'high');
    EEG.data=filtfilt(c,d,double(EEG.data));
    %plot(EEG.data(60000:120000)); hold on
    [b,a]=butter(2,30/(EEG.srate/2));
    EEG.data=filtfilt(b,a,double(EEG.data));
    %plot(EEG.data(60000:120000)); hold on
    
    tsec=(1:size(EEG.data,2))*(1/EEG.srate);
    filename=[save2folder '/filtecg' subid '_tms.txt'];
    fid=fopen(filename,'w');
    fprintf(fid, '%f %f \n', [tsec' EEG.data']');
    fclose(fid);
    clear EEG EEG
end
%% add fake cardiac phase events
load("/data/p_02186/TMS_ECG2/analyses/EEGsham/merge/cardiacphase_fake/fakeecg_equal_sysdys_sham.mat")
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_sham.set']);
    EEG= pop_select(EEG, 'nochannel', {'ECG','FDI'});
    
    allecg_sys=(sys_dys.faketrig_lat(sys_dys.subject==s &  sys_dys.systole==1)* EEG.srate)';
    allecg_dys=(sys_dys.faketrig_lat(sys_dys.subject==s &  sys_dys.diastole==1)* EEG.srate)';
    
    for i=1:length(allecg_sys)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='fake_sys';
        EEG.event(n_events+1).latency=(allecg_sys(i));
    end
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    for i=1:length(allecg_dys)
        n_events=length(EEG.event);
        EEG.event(n_events+1).type='fake_dys';
        EEG.event(n_events+1).latency=(allecg_dys(i));
    end
    EEG=eeg_checkset(EEG,'eventconsistency');
    EEG = pop_epoch(EEG, {'TMS'} , [-1.4 1]);
    EEG = pop_resample( EEG, 500);
    event_phase=EEG.event;
    
    save([save2folder '/' 'EEG_event_phase' subid], 'event_phase')
end