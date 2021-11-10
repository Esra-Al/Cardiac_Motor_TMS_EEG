clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
%rmpath('/home/raid1/esraal/Documents/eeglab2019_0/plugins/Fieldtrip-lite20200521')
eeglab
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGrest/';
save2folder = [eegpath 'merge']; mkdir(save2folder); cd(save2folder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
restcond={'rest'};
bl=1; %only one rest block
%%
for s = 1:subnum
    subid = subj_names{s};
    ALLEEG=[];
    EEG = pop_loadset([eegpath subj_names{s} '_' restcond{bl} '.set']);
    % EEG = pop_resample(EEG,fs); %don't resample
    EEG.chanlocs(31).labels='M2';
    EEG=pop_chanedit(EEG, 'lookup','/data/hu_esraal/Documents/eeglab2019_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc','eval','chans = pop_chancenter( chans, [],[]);');
 
    EEG.setname = subid;
    %create a new folder inside eegdata
    mkdir(subid); subfolder=[save2folder '/' subid];
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_rest.set'], 'filepath', subfolder);
end