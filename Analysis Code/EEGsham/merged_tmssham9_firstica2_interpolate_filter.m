clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
addpath('/data/hu_esraal/Documents/TMSEEG-4.0')
addpath('/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/')
eeglab
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGsham/';
eegfolder = [eegpath 'merge'];  cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
s=1
removeFrom  = -2;
removeTo = 15;
%%
for s = 36:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_sham_1_2_3_4_5.set']);

    %apply ica2 without filtering
    %copy ica weights
    EEGica=pop_loadset([subfolder subid '_tmsall_sham_1_2_3_4_5_6_7.set']);
    EEG.icasphere=EEGica.icasphere;
    EEG.icaweights=EEGica.icaweights;
    clear EEGica
    
    EEGcomps=pop_loadset([subfolder subid '_tmsall_sham_1_2_3_4_5_6_7_8.set']);
    EEG.comprej=EEGcomps.comprej;
    clear EEGcomps
    
    EEG=eeg_checkset(EEG,'ica');
    EEG = pop_subcomp(EEG,EEG.comprej, 0);
    EEG = eeg_checkset(EEG);
    
    EEG = tmseeg_addTMSTimeBack_interp(EEG, EEG.epoch_length);
    
    temp_removeFrom = find(EEG.times<=removeFrom);
    temp_removeTo = find(EEG.times<=removeTo);
    period2remove = temp_removeFrom(end)+1:temp_removeTo(end)-1;
    EEG = eeg_checkset(EEG);
    
    EEG.data(:,period2remove,:) = [];
    EEG.times(period2remove)    = [];
    EEG.TMS_period2remove_1     = period2remove;
    EEG.custom.times            = EEG.times;
    EEG.pnts                    = numel(EEG.times);
    EEG                         = eeg_checkset( EEG ); % some adjustments would be made
    nevents = length(EEG.event);
    for index = 1 : nevents
        EEG.event(index).latency=EEG.event(index).latency-(index-1)*length(period2remove);
    end
    
    EEG = tmseeg_addTMSTimeBack_interp(EEG, EEG.epoch_length);
      %Filter Design
    Fs=EEG.srate;
    ord = 2;
    [z1, p1, k1]  = butter(ord,[0.5 45]/(Fs/2),'bandpass');
    [xall1,yall2] = zp2sos(z1,p1,k1);
    notch_center=50; notch_size=5;
    [z2, p2, k2]  = butter(ord, [notch_center-(notch_size/2) notch_center+(notch_size/2)]/(Fs/2), 'stop');
    [xs1,xs2]     = zp2sos(z2,p2,k2); % Convert to 2nd order sections form
    
    %Apply Filter
    for ch=1:size(EEG.data,1)
        tempA=filtfilt(xall1,yall2,reshape(double(EEG.data(ch,:)),size(EEG.data,2),size(EEG.data,3)));
        tempB=filtfilt(xs1,xs2,double(tempA)); % apply notch filter
        EEG.data(ch,:,:)= double(tempB);
    end
    EEG = eeg_checkset( EEG );
    
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_firstica2_interpolate_filter.set'], 'filepath', subfolder, 'version', '7.3');
end