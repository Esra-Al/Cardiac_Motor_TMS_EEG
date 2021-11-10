function [ISI, jitter] = TMS_sequence_2stim(ISI, n_trials)
% 2 stimuli alternatingly; ISI in seconds
%  this function was used as TMS_sequence_2stim(2, 104) in our study
WaitSecs(25)
output_port = 23552; %LPT3 = 19456; LPT1(broken) = 888; LPT2 = 23552
trigger_length = 0.01; % 10 ms (maybe change this?)
key_idx = 27; % key index of ESC key

stim_sequence = repmat([17,18], 1, n_trials/2); % stimulation sequence starting with median stimulation
jitter = rand(n_trials, 1)-0.5;

for ii = 1:n_trials

    WaitSecs(ISI+jitter(ii)-trigger_length);
    lptwrite(output_port, stim_sequence(ii));
    WaitSecs(trigger_length);
    lptwrite(output_port, 0);

    % press ESC to stop script
    [~, ~, keyCode, ~] = KbCheck();
    if keyCode(key_idx) == 1 % press any key to stop script
        disp('Key pressed');
        break;
    end
    disp(ii)
end

end
