%% click calibration

% 1kHz tone

fs = 48000;
cf = 1000;
t = 0:1/fs:60-1/fs;
orig_gainf = 0.5;
ref_stim = orig_gainf*sin(2*pi*cf*t);

%plot(t,ref_stim)

%% click
click_dur = 0.0001; % 100 mu s
polarity = -1; %alternating polarity
rate = 9.1; % stim presentation rate  (9.1 Hz or 40.1)
stimL = 1/rate; % stim length (s)
nreps = 300; % repetitions (2000, 1000 in each polarity?)

RETSPL = 35.5;
peSPL = 80+RETSPL;


%gen stim
t = 0:1/fs:stimL-1/fs;
%click
polarity = -1;
click =polarity*[ones(round(click_dur*fs),1); ...
    zeros(round((stimL-click_dur)*fs)-1,1)];
click_stim = click;
click_stim = repmat(click,nreps,1);

%% play and adjust to wanted value (115.5dB)
ear = 2; % left or right
ppSPL = 115.5;
if ear ==2;
    gain = 10^((ppSPL-115.5)/20); % right
else
    gain = 10^((ppSPL-116)/20); % left
end
%gain = 1;
ref_stim = orig_gainf*sin(2*pi*cf*t);
ref_stim = ref_stim*gain;

%
%play 1kHz tone
KbName('UnifyKeyNames');
escKey = KbName('e');
quitKey = KbName('q');
lKey = KbName('LeftArrow');
try
    PsychPortAudio('GetOpenDeviceCount')
    PsychPortAudio('close');
catch
end
InitializePsychSound;
dev = PsychPortAudio('GetDevices');

%FOR PHY2
devid = 42; % 0 for PHYS2 HEAAUD
selectchannel =  [4 5;0 0]; %  [4 12;0 0]; ER2 + adat3
%devid = -1;
%selectchannel = [0 1;0 0];
nchans =  size(selectchannel,2);
pah = PsychPortAudio('Open', devid, [], 0, fs, nchans, [], [], selectchannel);
pa_status = PsychPortAudio('GetStatus',pah);
deviceId = pa_status.OutDeviceIndex;
PsychPortAudio('Volume',pah,0);
%PsychPortAudio('FillBuffer', pah, [stim_y trig_y]'); % chan x time
PsychPortAudio('FillBuffer', pah, [ref_stim' zeros(size(ref_stim'))]'); % chan x time
PsychPortAudio('Start', pah, 1, 0, 0, GetSecs+.1);
PsychPortAudio('Stop', pah, 1);
PsychPortAudio('Volume',pah,1);


try
    ListenChar(2)
    % get stim from stim struct
    %stim_y = stim.abr_stim{k};
    %trig_y = stim.trigger{k};
    
    % set trigger value
    %         if trig.is_connected()
    %            trig.set_trigger(stim.trigval(k));
    %         end
    
    %fill the buffer
    %PsychPortAudio('FillBuffer', pah,[stim_y trig_y]')
    if ear ==2
        PsychPortAudio('FillBuffer', pah,[zeros(size(ref_stim')) ref_stim']');
    else
        PsychPortAudio('FillBuffer', pah,[ref_stim' zeros(size(ref_stim'))]');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing screen
    clc
    disp(['Playing reference 1kHz tone'])
    
    %Get times and start PsychPortAudio
    offset = 0.2*rand;
    t_audio_start = GetSecs;
    PsychPortAudio('Start', pah, 1, t_audio_start+offset);
    % should not be started yet here due to offset
    s = PsychPortAudio('GetStatus', pah);
    t_audio = GetSecs;
    if ~s.Active
        while 1 % wait untill audio start
            s = PsychPortAudio('GetStatus', pah);
            [~,~,keyCode] = KbCheck(-1); % read key input
            if s.Active
                t_audio = GetSecs;
                break;
            end
            if keyCode(escKey)
                %fprintf('escape/n')
                PsychPortAudio('Stop', pah, 0);
                break;
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Audio start
    
    if s.Active % audio active
        while 1
            s = PsychPortAudio('GetStatus', pah);
            
            if ~s.Active
                break;
            end
            [~,~,keyCode] = KbCheck(-1); % read key input
            if keyCode
                fprintf('escape')
                PsychPortAudio('Stop', pah, 0);
                break;
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Sequence and response section
            while GetSecs<t_audio+length(ref_stim)/fs
                [keyIsDown,tResp,keyCode] = KbCheck;
                if keyIsDown
                    resp = find(keyCode); resp = resp(1); % responded key
                    
                    if resp==escKey
                        fprintf('escape')
                        PsychPortAudio('Stop', pah, 0);
                        break
                        
                    end
                    while KbCheck(-1);end
                    
                    
                elseif ~keyIsDown
                end
                
            end
        end
    end
    
    
    PsychPortAudio('Stop', pah, 1);
    clc
    while KbCheck; end
    WaitSecs(.1);
    KbWait(-1);
    [~,~,keyCode] = KbCheck;
    if find(keyCode)==quitKey
        commandwindow
        fprintf('escape')
        error('Exp end')
        
    else
    end
catch
    ListenChar(0)
    PsychPortAudio('Stop', pah);
    %     trig.disconnect();
    while KbCheck; end
    psychrethrow(psychlasterror);
end
ListenChar(0)
PsychPortAudio('Stop', pah);

clc
disp(['note the ppV value'])

%% play click and adjust to pp

%adjust to noted ppV
orig_gf = 0.5;
if ear ==2
    gain = 0.64; % right
else
    gain = 0.615; % left
end
click_stim_adjust = click_stim;
%gain = orig_gf*10^((peSPL-116.5)/20);
click_stim_adjust = gain*click_stim_adjust;


try
    PsychPortAudio('GetOpenDeviceCount')
    PsychPortAudio('close');
catch
end
InitializePsychSound;
dev = PsychPortAudio('GetDevices');

%FOR PHY2
devid = 42; % 0 for PHYS2 HEAAUD
selectchannel =  [4 5;0 0]; %  [4 12;0 0]; ER2 + adat3
%devid = -1;
%selectchannel = [0 1;0 0];
nchans =  size(selectchannel,2);
pah = PsychPortAudio('Open', devid, [], 0, fs, nchans, [], [], selectchannel);
pa_status = PsychPortAudio('GetStatus',pah);
deviceId = pa_status.OutDeviceIndex;
PsychPortAudio('Volume',pah,0);
%PsychPortAudio('FillBuffer', pah, [stim_y trig_y]'); % chan x time
if ear==2
    PsychPortAudio('FillBuffer', pah, [zeros(size(click_stim_adjust)) click_stim_adjust]'); % chan x time
else
    PsychPortAudio('FillBuffer', pah, [click_stim_adjust zeros(size(click_stim_adjust))]'); % chan x time
end
PsychPortAudio('Start', pah, 1, 0, 0, GetSecs+.1);
PsychPortAudio('Stop', pah, 1);
PsychPortAudio('Volume',pah,1);


try
    ListenChar(2)
    % get stim from stim struct
    %stim_y = stim.abr_stim{k};
    %trig_y = stim.trigger{k};
    
    % set trigger value
    %         if trig.is_connected()
    %            trig.set_trigger(stim.trigval(k));
    %         end
    
    %fill the buffer
    %PsychPortAudio('FillBuffer', pah,[stim_y trig_y]')
    PsychPortAudio('FillBuffer', pah,[click_stim_adjust click_stim_adjust]')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing screen
    clc
    disp(['Playing reference click'])
    
    %Get times and start PsychPortAudio
    offset = 0.2*rand;
    t_audio_start = GetSecs;
    PsychPortAudio('Start', pah, 1, t_audio_start+offset);
    % should not be started yet here due to offset
    s = PsychPortAudio('GetStatus', pah);
    t_audio = GetSecs;
    if ~s.Active
        while 1 % wait untill audio start
            s = PsychPortAudio('GetStatus', pah);
            [~,~,keyCode] = KbCheck(-1); % read key input
            if s.Active
                t_audio = GetSecs;
                break;
            end
            if keyCode(escKey)
                %fprintf('escape/n')
                PsychPortAudio('Stop', pah, 0);
                break;
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Audio start
    
    if s.Active % audio active
        while 1
            s = PsychPortAudio('GetStatus', pah);
            
            if ~s.Active
                break;
            end
            [~,~,keyCode] = KbCheck(-1); % read key input
            if keyCode
                fprintf('escape')
                PsychPortAudio('Stop', pah, 0);
                break;
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Sequence and response section
            while GetSecs<t_audio+length(click_stim_adjust)/fs
                [keyIsDown,tResp,keyCode] = KbCheck;
                if keyIsDown
                    resp = find(keyCode); resp = resp(1); % responded key
                    
                    if resp==escKey
                        fprintf('escape')
                        PsychPortAudio('Stop', pah, 0);
                        break
                        
                    end
                    while KbCheck(-1);end
                    
                    
                elseif ~keyIsDown
                end
                
            end
        end
    end
    
    
    PsychPortAudio('Stop', pah, 1);
    clc
    while KbCheck; end
    WaitSecs(.1);
    KbWait(-1);
    [~,~,keyCode] = KbCheck;
    if find(keyCode)==quitKey
        commandwindow
        fprintf('escape')
        error('Exp end')
        
    else
    end
catch
    ListenChar(0)
    PsychPortAudio('Stop', pah);
    %     trig.disconnect();
    while KbCheck; end
    psychrethrow(psychlasterror);
end
ListenChar(0)
PsychPortAudio('Stop', pah);
%54