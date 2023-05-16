function [tf]= field_tf(cond ,pb, label)

for s=1:length(pb)

    file=load(['TFA_' pb{s} cond '.mat']);
    tf.label=label;
    tf.freq=file.freq;
    tf.time=file.times/1000;
    tf.dimord='subj_chan_freq_time';
    tf.powspctrm(s,:,:,:)=file.TFdat;
            
end
end