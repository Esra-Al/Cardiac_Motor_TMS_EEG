clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
eeglab
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGtms/'; 
eegfolder = [eegpath 'merge'];  cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_default.set']);
    EEG= pop_select(EEG, 'nochannel', {'ECG','FDI'});
    EEG = pop_epoch(EEG, {'TMS'} , [-1.4 1]);
    EEG = pop_rmbase(EEG, [-110 -10]); 
    EEG.chanloc_orig = EEG.chanlocs; 
    EEG.epoch_length = length(EEG.times)/EEG.srate;
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_default_1.set'], 'filepath', subfolder, 'version', '7.3');
end