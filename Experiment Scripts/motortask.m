%%% Motor Task 
% Authors: Esra Al & Tilman Stephani, September 2019

fclose all;
close all;
cd D:\USER\Esra
addpath('D:\USER\Tilman\Other_toolboxes\Palamedes\Palamedes');
rand('state', sum(100*clock));
%% parameters

% Duration
t_resting_state_minutes = 3; % minutes
t_resting_state = t_resting_state_minutes * 60; % in seconds
task_dur = 3; % after start instruction and before end instruction to prevent visual effects

% triggers
tr_start = 101;
tr_end = 102;
tr_green=103;
tr_red=104;

% LPT port
sr = 5000; % of EEG amplifier; in Hz
port2 = 19456; % LPT port ID 
lptwrite(port2,0);

% Set screen properties
screenNumber        = max(Screen('Screens'));
screenProperties    = Screen('Resolution', screenNumber);
res                 = [screenProperties.width screenProperties.height];
pixDepth            = screenProperties.pixelSize;
nbBuffers           = 2;
backColor           = [0 0 0]; % black
[window, winRect]   = Screen('OpenWindow', screenNumber, backColor, [], pixDepth, nbBuffers); % open window
[P.winCenterX, P.winCenterY] = WindowCenter(window); % center of the screen
HideCursor;

% Define text format
txtsize = 18; 
fixsize = 60;
txtspace = 2; % distance between lines
Screen('TextFont', window, 'Arial');
Screen('TextSize', window, txtsize);
Screen('TextStyle', window, 0); % letter width normal

% Color specifications
P.grey = [140 140 140];
P.white = [255 255 255];
P.green = [0 255 0];
P.red = [255 0 0];

% Fixation cue
P.cueRect = [0 0 40 40];
P.cuePos = CenterRectOnPoint(P.cueRect, P.winCenterX, P.winCenterY);
P.cueWidth = 15;
%% start instructions

instr = ['Willkommen! \n\n' ...
    'Dieser Teil enhält eine Motoraufgabe.\n' ...
    'Bitte entspannen Sie sich und schauen Sie auf den Fixationskreis \n' ...
    'in der Mitte vom Bildschirm.\n\n' ...
    ' Wenn der Kreis rot ist, drücken Sie den Drucksensor\n' ...
    'Wenn der Kreis grün ist, lassen Sie ihn los\n' ...
    'Diese Messung wird ' num2str(t_resting_state_minutes) ' Minuten dauern. \n\n'];
DrawFormattedText(window, instr, 'center', 'center', P.white, [], [], [], txtspace); % put text in hidden screen
Screen('flip', window); % show hidden screen
[secs, keyCode] = KbStrokeWait(); % wait for button press
if find(keyCode==1) == 27 % if ESC has been pressed
    ShowCursor;
    sca; % close presentation screen
end
%% Start
% end trigger
lptwrite(port2, tr_start); % send trigger
WaitSecs(1/sr*2); % make sure trigger is not missed
lptwrite(port2, 0);
WaitSecs(task_dur)
for c=1:60
    instr = 'o';
    if mod(c,2)==1
        Screen('TextSize', window, fixsize); % increase font size for fixation cross
        DrawFormattedText(window, instr, 'center', 'center', P.red, [], [], [], txtspace); % put text in hidden screen
        Screen('flip', window); % show hidden screen

        

        % start trigger
        lptwrite(port2, tr_red); % send trigger
        WaitSecs(1/sr*2); % make sure trigger is not missed
        lptwrite(port2, 0);
    else
        
        Screen('TextSize', window, fixsize); % increase font size for fixation cross
        DrawFormattedText(window, instr, 'center', 'center', P.green, [], [], [], txtspace); % put text in hidden screen
        Screen('flip', window); % show hidden screen


        % start trigger
        lptwrite(port2, tr_green); % send trigger
        WaitSecs(1/sr*2); % make sure trigger is not missed
        lptwrite(port2, 0);
    end
    % hold ESC for quitting
    [keyIsDown, secs, keyCode, deltaSecs]  = KbCheck(); % check for button press
    if keyIsDown
        if find(keyCode==1) == 27 % if ESC has been pressed
            break % stop sequence
        end
    end

   WaitSecs(task_dur);
end
% end trigger
lptwrite(port2, tr_end); % send trigger
WaitSecs(1/sr*2); % make sure trigger is not missed
lptwrite(port2, 0);

% "buffer" interval
WaitSecs(task_dur);
%% End screen
instr = ['Messung beendet! \n\n'...
    'Bitte Versuchsleiter Bescheid geben.'];
Screen('TextSize', window, txtsize); % decrease font size for text
DrawFormattedText(window, instr, 'center', 'center', P.white, [], [], [], txtspace); % put text in hidden screen
Screen('flip', window); % show hidden screen
[secs, keyCode] = KbStrokeWait(); % wait for button press
if find(keyCode==1) == 27 % if ESC has been pressed
    ShowCursor;
    sca; % close presentation screen
end

ShowCursor;
sca;






