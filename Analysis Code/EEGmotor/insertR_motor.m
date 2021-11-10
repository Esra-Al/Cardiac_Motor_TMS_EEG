clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
eeglab
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGmotor/';
cd(eegpath)
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
restcond={'motor'};
%trigger 104 - red dot- press and 103: green dot- release-relax
%%
EEG = eeg_checkset(EEG,'eventconsistency');
%pop_saveset(EEG, [eegpath subj_names{s} '_' restcond{bl} '.set'])
%%
s=36;
bl=1;
subid=subj_names{s};
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
clear EEG HEP
EEG=pop_loadset([eegpath subj_names{s} '_' restcond{bl} '.set']);
%
ecg_chn =32; %ECG electrode
ecg = double(EEG.data(ecg_chn,:));
ecg = ecg(:); % make sure ecg is one-column vector
[c,d]=butter(2,0.5/(EEG.srate/2),'high');
[b,a]=butter(2,30/(EEG.srate/2));
ecg=filtfilt(c,d,double(ecg));
ecg=filtfilt(b,a,double(ecg));
HEP.ecg = (ecg-min(ecg))/(max(ecg)-min(ecg));
HEP.srate = EEG.srate;
HEP.winsec = 22;
HEP.sec_ini = 0;
HEP.qrs = heplab_fastdetect(ecg,EEG.srate);
clear ecg ecg_chn
heplab
%%
EEG = eeg_checkset(EEG,'eventconsistency');
pop_saveset(EEG, [eegpath subj_names{s} '_' restcond{bl} '.set'])
%% check whether ECG event saved
ecg_detected = NaN(subnum, 1);
for s = 1:subnum
    subid = subj_names{s};
    %for bl=1:4
        EEG=pop_loadset([eegpath subj_names{s} '_' restcond{bl} '.set']);
        [EEGout,ind] = pop_selectevent(EEG,'type','ECG');
        if ~isempty(ind)
            ecg_detected(s,bl) = length(ind);
        end
        clear EEG
    %end
end
save ecg_events ecg_detected