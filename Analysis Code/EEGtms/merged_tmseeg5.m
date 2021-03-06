%remove and plot ICA1 components
%%
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
xmin = -150;
xmax = 150;

for s = 1:subnum
    clear EEG I
    subid = subj_names{s}; subfolder=[eegfolder '/' subid '/'];
    EEG = pop_loadset([subfolder subid '_tmsall_default_1_2_3_4.set']);
    times = EEG.times;
    % calculating the variance explained by the components
    mact = mean(eeg_getdatact(EEG,'component',1:size(EEG.icawinv,2)),3);
    twin = times > xmin & times < xmax;
    
    for k=1:size(EEG.icawinv,2)
        a = EEG.icawinv(:,k)*mact(k,:);   %project components
        v(k) = max(abs(max(a(:,twin),[],2)-min(a(:,twin),[],2))); %#ok %variance
    end
    [~, I] = sort(v,'descend'); % sorting in order of variance
    
    y_temp = mean(EEG.data,3);
    temp_times = min(times):1000/EEG.srate:max(times);
    seg = find(temp_times>xmin & temp_times<xmax);
    x = temp_times(seg);
    %Insert NaN values to fill space where TMS pulse was removed
    ix = min(EEG.TMS_period2remove_1);
    EEG_temp = y_temp;
    rm_pulse_fill = NaN(size(EEG_temp,1),length(EEG.TMS_period2remove_1));
    y_temp = cat(2,EEG_temp(:,1:ix-1),rm_pulse_fill,EEG_temp(:,ix:end));
    y_bef = y_temp(:,seg)';
    %sel=3; %remove first  3 components
    for sel=3:4 %remove first 3 or 4 components
        selected=1:sel;
        projection = eeg_getdatact(EEG,'rmcomps',I(selected));
        temp_times = min(times):1000/EEG.srate:max(times);
        
        proj_temp = mean(projection,3);
        %Insert NaN values to fill space where TMS pulse was removed
        proj_temp = cat(2,proj_temp(:,1:ix-1),rm_pulse_fill,proj_temp(:,ix:end));
        y = proj_temp(:,seg)';
        
        %Plot Data channels overlain, averaged across trials
        figure;
        subplot(1,2,1)
        plot(x,y_bef); xlim([xmin xmax])
        title('Before component removal')
        xlabel('Time (ms)')
        ylabel(['Amplitude (' char(0181) 'V)'])
        
        subplot(1,2,2)
        plot(x,y); xlim([xmin xmax])
        title('After component removal')
        xlabel( 'Time (ms)')
        ylabel(['Amplitude (' char(0181) 'V)'])
        savefig([subfolder subid 'bef_aft_ICA1_first' num2str(sel) 'comp.fig']) 
        close all
    end
    %Save, remove 3 first components
    selected=1:3;  %remove 3 first components
    EEG.decay_icawinv     = EEG.icawinv;
    EEG = pop_subcomp( EEG, I(selected), 0); % remove components
    EEG.decaycomp_removed = I(selected);
    disp(I(selected))
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG, 'filename', [subid '_tmsall_default_1_2_3_4_5.set'], 'filepath', subfolder, 'version', '7.3');
end