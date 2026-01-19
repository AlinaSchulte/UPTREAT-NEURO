function ffr_SW_main(dat)
clc
addpath(genpath(pwd))

%% script to run FFRs


%% Subid
subid = dat.subject;
%% Choose calibration file
%% Choose calibration file
ears = {'L','R'};
choiceTT.ear = ears(dat.ear);
choiceTT.ear_noise = ears(setxor(dat.ear,[1,2]));
choiceTT.filename = 'ER3_MAG_48k_current';

%% Generate tone stim
fs = 48000;
%fs_eeg = 64;

% init parameeters
%tonel = 326;
stim = struct;

stimL = 3;
stimrate = 3;
fmod = 2;

%conditions:
% reps x conds x Nperiods x 2s
reps = 1;
conds = [];
conds = repmat([1:2],1,reps);
%conds = Shuffle(conds);
Nperiods = 85; % periods per trial x 6 stims

% total length of experiment:
tl = reps*length(conds)*stimL*Nperiods/60;
fprintf('Total length of experiment %.2f min\n',tl);
%%
%tic
% generating trials for each condition
stim = {};
stim.fc = [repmat(326,1,2)];
stim.polarity = [repmat([1,-1],1)];
stim.noise = [repmat([0 0],1,1)];
stim.modf = repmat(4,1,2);
stim.ear = repmat(dat.ear,1,2);
stim.targetlevel = repmat(85,1,2);

for tt=1:length(conds)
    cc=conds(tt);     % condition
    trig_val = conds(tt)*10; % Trigger value
    
    tone_c = [];trig_c = [];
    % generate a trial
    
    %generate tone stimuli
    f_tone = stim.fc(cc);
    t=0:1/fs:stimL-1/fs;
    tones = sin(2*pi*f_tone*t);
    
     % Calibration right channel
    choiceTT.targetlevel = stim.targetlevel(cc)+nHLfilter(stim.targetlevel(cc),f_tone); %nHL compensation (dB)
    choiceTT.signal = tones;
    tones = calfilter(choiceTT.ear, choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs); 
    
    % add modulation
    tonemod = zeros(1,stimL*fs);
    ons = 1:fs/fmod:stimL*fs+fs/fmod;
    Lr = .01*fs; % ramp length
    r = [sin(linspace(0, pi/2, Lr)) ones(1,fs/fmod/2-Lr*2) sin(linspace(pi/2, 0, Lr))]'; % ramp

    for ii = 1:length(ons)-1
        tonemod(ons(ii):ons(ii)+fs/fmod/2-1) = r;
    end
    tones = tones'.*tonemod;% apply slow on/off modulation (6 tones)

    % trigger pulse

    %pulse every tone burst
    tpulse = [ones(1,round(1e-3 * fs)), zeros(1,ons(2)-1-round(1e-3 * fs))]';
    tpulse = repmat(tpulse,fmod*stimrate,1);
    
    
    % add silence (~1 s)
    sil_p = .6;
    tones_tmp = [tones  zeros(1,sil_p*fs)];
    tpulse = [tpulse;zeros(sil_p*fs,1)];
    %repmat out to full length
    tones = repmat(tones_tmp',Nperiods,1)*stim.polarity(cc);
    trig = repmat(tpulse,Nperiods,1);
    
    %apply 5s relax ramp to begining and end
    RAMP = ones(size(tones,1)-length(tones_tmp)*4,1);
    rin =([0:1/(length(tones_tmp)*2):1-1/(length(tones_tmp)*2)]).^2;
    RAMP = [rin';RAMP;flip(rin)'];
    
    %apply the ramp
    tones= RAMP.*tones;
    trig = [zeros(size(rin'));trig(1:end-2*length(rin));zeros(size(rin'))];
    % add noise?
    if stim.noise(cc)
        %TEN noise
        noiseL = length(tones)/fs;
        flow = 20;fhigh=20000;
        [noise,~,~] = TEN_noise(noiseL,fs,flow,fhigh);
        
        %adding noise
        pt =sum(abs(tones.^2)); %tone power
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
clc
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
PsychPortAudio('FillBuffer', pah, [stim_y trig_y]'); % chan x time
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
disp(['length of experiment: ',num2str(sum(ltmp)/60),' min'])
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
        PsychPortAudio('FillBuffer', pah,[stim_y trig_y]')
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
    save(['ffr_SW_stim_',subid,'_aborted.mat'],'stim','-v7.3');
    cd ..
    psychrethrow(psychlasterror);
end
ListenChar(0);
trig.disconnect();
delete(trig);
cd('./_data')
stim.tones = [];
save(['ffr_SW_stim_',subid,'.mat'],'stim','-v7.3');
cd ..
disp('Experiment Done!')
%
end
