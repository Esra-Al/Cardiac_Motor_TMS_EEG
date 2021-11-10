
%run ICA2 components
%%
clc; clear;
addpath('/data/hu_esraal/Documents/FastICA_25')
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
eeglab
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGsham/';
eegfolder = [eegpath 'merge'];  cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);

for s = 1:subnum
    clear EEG
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsalsham_1_2_3_4_5_6.set']);
    h1 = msgbox('Running ICA2,now!');
    EEG = pop_runica(EEG, 'icatype' ,'fastica','g','tanh','approach','symm');
    EEG = eeg_checkset(EEG);
    close(h1);
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_sham_1_2_3_4_5_6_7.set'], 'filepath', subfolder, 'version', '7.3');
end
%%
% addpath('/data/p_02186/TMS_ECG2/analyses/TESA1.1.1/')
% for s = 1:subnum
%     clear EEG
%     subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
%     EEG = pop_loadset([subfolder subid '_tmsall_sham_1_2_3_4_5_6_7.set']);
%     EEG_tes = pop_tesa_compselect( EEG,'compCheck', 'off','comps',[],'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','on','muscleThresh',0.6,'muscleFreqIn',[30 100],'muscleFeedback','off','elecNoise','off','elecNoiseThresh',6,'elecNoiseFeedback','off' );
%     EEG.ica2_rejectTmsMuscle = EEG_tes.icaCompClass.TESA2.rejectTmsMuscle;
%     EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_sham_1_2_3_4_5_6_7.set'], 'filepath', subfolder, 'version', '7.3');
% end