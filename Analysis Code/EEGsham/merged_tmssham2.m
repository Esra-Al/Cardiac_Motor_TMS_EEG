clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
eeglab
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGsham/';
eegfolder = [eegpath 'merge'];  cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
removeFrom  = -2;
removeTo = 8;

for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_sham_1.set']);
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
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_sham_1_2.set'], 'filepath', subfolder, 'version', '7.3');
end