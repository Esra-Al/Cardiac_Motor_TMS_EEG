%remove and plot ICA1 components
%%
clc; clear;
addpath('/data/hu_esraal/Documents/eeglab2019_0/')
eeglab
%%
eegfolder= '/data/p_02186/TMS_ECG2/analyses/EEGsham/merge'; cd(eegfolder);
subj_names = {'VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11', ...
    'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23', ...
    'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37'};
subnum=length(subj_names);
xmin = -150;
xmax = 150;
removeFrom  = -2;
removeTo = 8;
%%
for s = 1:subnum
    clear EEG I
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_sham_1_2_3_4.set']);
    EEGorig=EEG;
%     times = EEG.times;
%     % calculating the variance explained by the components
%     mact = mean(eeg_getdatact(EEG,'component',1:size(EEG.icawinv,2)),3);
%     twin = times > xmin & times < xmax;
%     
%     for k=1:size(EEG.icawinv,2)
%         a = EEG.icawinv(:,k)*mact(k,:);   %project components
%         v(k) = max(abs(max(a(:,twin),[],2)-min(a(:,twin),[],2))); %#ok %variance
%     end
%     [~, I] = sort(v,'descend'); % sorting in order of variance
%     
%     y_temp = mean(EEG.data,3);
%     temp_times = min(times):1000/EEG.srate:max(times);
%     seg = find(temp_times>xmin & temp_times<xmax);
%     x = temp_times(seg);
%     
%     selected=1:3;  %remove 3 first components
%     comprej= I(selected);

    ECG = pop_loadset([subfolder subid '_tmsall_sham.set']);
    ECG= pop_select(ECG, 'channel', {'ECG'});
    ECG = pop_epoch(ECG, {'TMS'} , [-1.4 1]);
    ECG = pop_rmbase(ECG, [-110 -10]);
    
    temp_removeFrom = find(ECG.times<=removeFrom);
    temp_removeTo = find(ECG.times<=removeTo);
    period2remove = temp_removeFrom(end)+1:temp_removeTo(end)-1;
    ECG = eeg_checkset( ECG );
    ECG.data(:,period2remove,:) = [];
    ECG.times(period2remove)    = [];

    EEG.data(end+1,:,:) = ECG.data;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'ECG';
    EEG = eeg_checkset(EEG);
    
    ncomp= size(EEG.icawinv,2); % ncomp is number of components
    icaacts = eeg_getdatact(EEG,'component',1:ncomp);
    EEG.icaact = icaacts;
    EpochSecs = 2.4;
    EEG3 = EEG; EEG3.data=EEG3.icaact;
    
    %
    chanEEG = EEG.data(64,:)';
    ICs = EEG3.data(:,:)'; 
    c  = abs(corr(ICs,chanEEG))';
    corthresh = mean(c)+std(c)*1.5;
    rej = c > corthresh ; 
    [B,ecgrej] = maxk(c,2);
    
    %comprej=[comprej ecgrej];
    
    EEG = pop_subcomp(EEGorig, ecgrej, 0); % remove components
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_sham_1_2_3_4_5_ecg.set'], 'filepath', subfolder, 'version', '7.3');
end