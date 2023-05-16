function [tf]= field_tf(cond ,pb,subnum, label)

for s=1:subnum

    file=load(['TFA_VP' pb(s).name cond '.mat']);
    tf.label=label;
    tf.freq=file.freq;
    tf.time=file.times/1000;
    tf.dimord='subj_chan_freq_time';
    tf.powspctrm(s,:,:,:)=file.TFdat;
            
end
end