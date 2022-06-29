function [EEG] = tmseeg_addTMSTimeBack_interp(EEGkk, EpochSecs)

EEG = EEGkk;
discon_low = [];
discon_high = [];

for i=1:(size(EEGkk.times,2)-1) % find discontinuity in times array
    
    if( (EEGkk.times(i+1)-(EEGkk.times(i)+1000/EEGkk.srate)) >= (1000/EEGkk.srate))
        discon_low = [discon_low; EEGkk.times(i)]; %#ok
        discon_high = [discon_high; EEGkk.times(i+1)]; %#ok
    end
    
end

discon = [discon_low discon_high];

for i=1:size(discon,1) % sometimes finds a discontinuity at end of epoch (i.e., 998ms -> 999.0000ms for epochs (-1000ms to +1000ms))
    
    if(discon(i,1)>EEGkk.times(size(EEGkk.times,2)-5))
        discon(i,:) = [];
    end
    
end

ch = size(EEGkk.data,1);
tr = size(EEGkk.data,3);
EEG.data = zeros(ch,EEGkk.srate*EpochSecs,tr);


sp_low = discon(1,1);
sp_high = discon(1,2);


% single pulse paradigm
w1 = find(EEGkk.times==sp_low);
w2 = find(EEGkk.times==sp_high);
gap = ceil((sp_high - sp_low)*(EEGkk.srate/1000) - 1);

EEG.pnts = EEGkk.pnts + gap ;
EEG.data(:,1:w1,:) = EEGkk.data(:,1:w1,:);
EEG.data(:,w1+gap+1:EEG.pnts,:) = EEGkk.data(:,w2:EEGkk.pnts,:);
EEG.data(:,w1+1:w1+gap,:) = NaN;
EEG.times = floor(EEGkk.times(1)):(1000/EEGkk.srate):EEGkk.times(end);

nevents = length(EEG.event);
for index = 1 : nevents
    EEG.event(index).latency=EEG.event(index).latency+(index-1)*gap;
end

%interpolate
prewindow  = round(0.01 * EEG.srate);  % Express window in samples
postwindow = round(0.01* EEG.srate); % Express window in samples

tim=EEG.times;
ntrl=EEG.trials;
nchan=EEG.nbchan;
for tr=1:ntrl
    for ch=1:nchan
        dat=EEG.data(ch,:,tr);
        replace = isnan(dat); % Find samples that have been replaced by nans
        onset  = find(diff([0 replace])>0);
        offset = find(diff([replace 0])<0);
        for k=1:numel(onset)
            begsample = onset(k)-prewindow;
            endsample = offset(k)+postwindow;
           
            x = tim(begsample:endsample);
            y = dat(begsample:endsample);
            xx = x; % this is where we want to know the interpolated values
            x = x(~replace(begsample:endsample)); % remove the part that needs to be interpolated
            y = y(~replace(begsample:endsample)); % remove the part that needs to be interpolated
            
             yy = interp1(x, y, xx, 'pchip'); % this may contain nans
            
            % The default extrapolation behavior of INTERP1 with four input arguments is to
            % extrapolate for 'spline', 'pchip' and 'makima', and to use nan for other
            % methods.
            
            if begsample==1
                % there may be nans at the beginning, replace the data with mean of the values that are not nan
                f = find(~isnan(yy), 1, 'first');
                yy(1:f-1) = nanmean(yy);
            elseif endsample==numel(dat)
                % there may be nans at the end, replace the data with mean of the values that are not nan
                f = find(~isnan(yy), 1, 'last');
                yy(f+1:end) = nanmean(yy);
            end
            
            % insert the interpolated data
            EEG.data(ch,begsample:endsample,tr) = yy;
            
        end % for all nan-segments
        
    end
end

