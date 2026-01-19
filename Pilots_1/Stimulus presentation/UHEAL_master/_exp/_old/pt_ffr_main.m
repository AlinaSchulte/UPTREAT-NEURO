function pt_ffr_main(dat)
clc
addpath(genpath(pwd))

%% script to run FFRs
%1) 4Hz 706Hz quiet
%2) 4Hz 706Hz noise contra (5dB SNR)
%3) 4Hz 242Hz quiet
%4) 4Hz 242Hz noise contra (5dB SNR)
%5-8) 1-4 in alternating polarities

%8 conditions in total

%% Subid
subid = dat.subject;
%% Choose calibration file
ears = {'L','R'};
choiceTT.ear = ears(dat.ear);
choiceTT.ear_noise = ears(setxor(dat.ear,[1,2]));
choiceTT.filename = 'ER3_MAG_48k_current';

%% Generate tone stim
fs = 48000;
%fs_eeg = 64;

% init parameeters
tonel = .5*fs;
stim = struct;

%conditions:
% reps x conds x Nperiods x 2s
reps = 1;
conds = [];
conds = repmat([1:8],1,reps);
conds = Shuffle(conds);
Nperiods = 500; % periods per trial

% total length of experiment:
tl = reps*length(conds)*tonel/fs*Nperiods/60;
fprintf('Total length of experiment %.2f min\n',tl);
%%
%tic
% generating trials for each condition
stim = {};
stim.fc = [repmat(706,1,4) repmat(326,1,4)];
stim.polarity = [repmat([1,-1],1,4)];
stim.noise = [repmat([0 0 1 1],1,2)];
stim.modf = repmat(4,1,8);
stim.ear = repmat(dat.ear,1,8);
stim.targetlevel = repmat(85,1,8);

for tt=1:length(conds)
    cc=conds(tt);     % condition
    trig_val = conds(tt)*10; % Trigger value

    
    tone_c = [];trig_c = [];
    % generate a trial
    
    %generate tone stimuli
    f_tone = stim.fc(cc);
    t=0:1/fs:tonel/fs-1/fs;
    tones = sin(2*pi*f_tone*t);
    
    % add modulation
    moddepth = 1;
    modl = stim.modf(cc);
    ttt = 0:1/fs:length(tones)/fs-1/fs;
    mod= (1+ (moddepth*(sin(2*pi*modl*ttt - pi/2))));
    mod = mod/max(abs(mod));
    tones = mod.*tones;
    tones = tones';

    % trigger pulse
    tpulse = [ones(1,round(1e-3 * fs)) zeros(1,tonel-round(1e-3 * fs))]';
    
    tones = repmat(tones,Nperiods,1)*stim.polarity(cc);
    trig = repmat(tpulse,Nperiods,1)';
    
    % Calibration right channel
    choiceTT.targetlevel = stim.targetlevel(cc)+nHLfilter(stim.targetlevel(cc),f_tone); %nHL compensation (dB)
    choiceTT.signal = tones;
    tones = calfilter(choiceTT.ear, choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
    % cal for noise ear
    tones_n = calfilter(choiceTT.ear_noise, choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
    
    
    % add noise?
    if stim.noise(cc)
        %TEN noise
        noiseL = length(tones)/fs;
        flow = 20;fhigh=20000;
        [noise,~,~] = TEN_noise(noiseL,fs,flow,fhigh);
        
        %adding noise
        pt =sum(abs(tones_n.^2)); %tone power
        %5dB SNR noise
        SNRc = 5;
        SNRtmp=pt/(10^(SNRc/10)); %required noise power
        pnoise = sum(abs(noise.^2)); %noise power
        tmpnoise = (noise'./sqrt(pnoise)).*sqrt(SNRtmp);
    else
        tmpnoise = zeros(size(tones));
    end
    
    
    % Which ear
    if strcmp(choiceTT.ear, 'R')
        stim.tones{tt} = [tmpnoise tones];
    else
        stim.tones{tt} = [tones tmpnoise];
    end
    
    % Make struct
    stim.trigger{tt} = trig;
    stim.time = [0:1/fs:length(tones)/fs-1/fs]';
    stim.trigval(tt) = trig_val;
    stim.conds = conds;
    
    
    clc
    disp(num2str(tt))
    
end
%toc
%reshuffling
stim.fc = stim.fc(conds);
stim.polarity = stim.polarity(conds);
stim.noise = stim.noise(conds);
stim.modf = stim.modf(conds);
stim.ear = stim.ear(conds);
stim.targetlevel = stim.targetlevel(conds);
stim.trigval = stim.trigval';
stim.experiment_time = dat.experiment_time;
stim.dat = dat;

%%
ltmp = [];
for ii = 1:length(stim.tones)
    ltmp(ii) = length(stim.tones{ii})/fs;
    %soundsc(stim.tones{ii},fs)
    %pause()
end
disp(['length of experiment: ',num2str(sum(ltmp)/60),' min'])

%% Experimental part
stim_y = stim.tones{1};
%add trigger here
trig_y = stim.trigger{1};
try
    PsychPortAudio('GetOpenDeviceCount')
    PsychPortAudio('close');
catch
end
InitializePsychSound;
dev = PsychPortAudio('GetDevices');

%FOR PHY2
devid = 42; % 0 for PHYS2 HEAAUD
selectchannel =  [4 5 17;0 0 0]; %  [4 12;0 0]; ER2 + adat3
%devid = -1;
%selectchannel = [0 1;0 0];
nchans =  size(selectchannel,2);
pah = PsychPortAudio('Open', devid, [], 0, fs, nchans, [], [], selectchannel);
pa_status = PsychPortAudio('GetStatus',pah);
deviceId = pa_status.OutDeviceIndex;
PsychPortAudio('Volume',pah,0);
PsychPortAudio('FillBuffer', pah, [stim_y trig_y']'); % chan x time
%PsychPortAudio('FillBuffer', pah, [stim_y']); % chan x time
PsychPortAudio('Start', pah, 1, 0, 0, GetSecs+.1);
PsychPortAudio('Stop', pah, 1);
PsychPortAudio('Volume',pah,1);

% init triggerbox
trig = HEATriggerbox();
trig.find_triggerbox_win();
trig.connect();
if trig.is_connected()
   trig.set_trigger(stim.trigval(1));
end

%% init kB
KbName('UnifyKeyNames');
escKey = KbName('e');
quitKey = KbName('q');

%% Experimental loop
clc
disp('Ready')
while KbCheck; end
WaitSecs(.1);
KbWait(-1);

Ntrials = length(stim.tones);
resp_mat=[];
try
    for k=1:Ntrials
        ListenChar(2)
        % get stim from stim struct
        stim_y = stim.tones{k};
        trig_y = stim.trigger{k};
        
        % set trigger value
                if trig.is_connected()
                   trig.set_trigger(stim.trigval(k));
                end
        
        %fill the buffer
        PsychPortAudio('FillBuffer', pah,[stim_y trig_y']')
        %PsychPortAudio('FillBuffer', pah,stim_y')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Drawing screen
        clc
        disp(['Playing stimulus '  num2str(stim.fc(k)) ' Hz. Noise ' num2str(stim.noise(k))  ', Polarity ' num2str(stim.polarity(k))])
        
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
                while GetSecs<t_audio+length(stim_y)/fs
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
        disp(['Trial ' num2str(...,
            [num2str(k), '/' ...,
            ,num2str(Ntrials)])])
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
        
        
    end
catch
    
    PsychPortAudio('Stop', pah);
    ListenChar(0);
    trig.disconnect();
    delete(trig);
    while KbCheck; end
    
    cd('./_data')
    stim.tones = [];
    save(['pt_ffr_stim_',subid,'_aborted.mat'],'stim','-v7.3');
    cd ..
    psychrethrow(psychlasterror);
end
ListenChar(0);
trig.disconnect();
delete(trig);
cd('./_data')
stim.tones = [];
save(['pt_ffr_stim_',subid,'.mat'],'stim','-v7.3');
cd ..
disp('Experiment Done!')
%
end
