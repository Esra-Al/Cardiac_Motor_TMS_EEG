EpochSecs = 2.4;
EEG2 = tmseeg_addTMSTimeBack(EEG, EpochSecs);
EEG3 = EEG; EEG3.data=EEG3.icaact;

%
ncomp= size(EEG.icawinv,2); % ncomp is number of components
chanEEG = EEG.data(64,:)';
ICs = EEG3.data(:,:)'; 
c  = abs(corr(ICs,chanEEG))';
corthresh = mean(c)+std(c)*1.5;
rej = c > corthresh ; 
[B,rej_comps] = maxk(c,3) ;
%

EEG3 = tmseeg_addTMSTimeBack(EEG3, EpochSecs);
ecginds=find(strcmp('ECG',{EEG2.event.type}));
ecglat=[EEG2.event(ecginds).latency];
tmslat=[EEG2.event(strcmp('TMS',{EEG2.event.type})).latency];

for st=1:length(tmslat)
    pre_ecg= find(tmslat(st)>ecglat, 1, 'last');
    rind1=ecginds(pre_ecg);
    
    if (tmslat(st)-ecglat(pre_ecg))/EEG2.srate >0.4 && (tmslat(st)-ecglat(pre_ecg))/EEG2.srate <1.4
        EEG2.event(rind1).type= 'pre_ECG';
    end
end

EEG2_ecg = pop_select(EEG2, 'channel', {'ECG'});
EEG2_ecg = pop_epoch(EEG2_ecg, {'pre_ECG'} , [-0.05 0.4]);
EEG2_ica=EEG2; 
EEG2_ica.data=EEG3.data;
EEG2_ica = pop_epoch(EEG2_ica, {'pre_ECG'} , [-0.05 0.4]);
%%
disp(rej_comps)
for r=1:length(rej_comps)
    set(groot,'defaultLineLineWidth',3.0) 
    comp = rej_comps(r);
    figure
    plot(EEG2_ica.times, mean(EEG2_ica.data(comp,:,:),3))
    hold on
    plot(EEG2_ecg.times, mean(EEG2_ecg.data(1,:,:),3)/300)
    legend('IC', 'ecg', 'FontSize',15)
    title(num2str(comp))
end


% disp(['Muscle components ###########  ' num2str(EEG.ica2_rejectTmsMuscle) '   ###########'])
