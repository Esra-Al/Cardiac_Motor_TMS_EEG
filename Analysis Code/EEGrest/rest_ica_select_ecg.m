ncomp= size(EEG.icawinv,2); % ncomp is number of components
chanEEG = EEG.data(64,:)';
ICs = EEG.icaact(:,:)'; 
c  = abs(corr(ICs,chanEEG))';
corthresh = mean(c)+std(c)*1.5;
rej = c > corthresh ; 
[B,rej_comps] = maxk(c,3) ;
%
%%
disp(rej_comps)
for r=1:length(rej_comps)
    set(groot,'defaultLineLineWidth',3.0) 
    comp = rej_comps(r);
    figure
    plot(EEG.times, mean(EEG.icaact(comp,:,:),3))
    hold on
    plot(EEG.times, mean(EEG.data(64,:,:),3)/1000)
    legend('IC', 'ecg', 'FontSize',15)
    title(num2str(comp))
end