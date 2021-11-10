clc; clear;
addpath('/data/hu_esraal/Documents/FastICA_25')
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
addpath('/data/pt_02333/esra_codes/')
eeglab
%%
eegpath= '/data/p_02186/TMS_ECG2/analyses/EEGrest/';
eegfolder = [eegpath 'merge'];  cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
grandtrials=cell(subnum,1); epochnum=NaN(subnum,1);
cfg2=[]; cfg2.baseline=[-0.1 0];

%%
for s = 1:subnum
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG =  pop_loadset([subfolder subid '_tmsall_rest.set']);
    EEG= pop_select(EEG, 'nochannel', {'ECG','FDI'});
    
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
    
    %h1 = msgbox('Running ICA,now!');
    EEG = pop_runica(EEG, 'icatype' ,'fastica','g','tanh','approach','symm');
    EEG = eeg_checkset(EEG);
    %close(h1);
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_rest_ica.set'], 'filepath', subfolder, 'version', '7.3');
end
