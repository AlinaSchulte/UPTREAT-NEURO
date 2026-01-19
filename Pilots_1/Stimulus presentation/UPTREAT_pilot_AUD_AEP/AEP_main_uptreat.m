% could not find main script for previous function. made main script out of
% this
clc
addpath(genpath(pwd))
 
%% script to run AEP
 
%8 conditions in total

    % load stims or subject files?? where is dat and choiceTT coming from?
%% Subid
subid = dat.subject;
%% Choose calibration file
ears = {'L','R'};
choiceTT.ear = ears(dat.ear);
choiceTT.ear_noise = ears(setxor(dat.ear,[1,2]));
choiceTT.filename = 'ER3_MAG_48k_current';
 
%% Generate tone stim
% stimulus parameters # adjusted to uptreat
%% General stimulus parameters
stim = struct;
stim.fs = 48000;
stim.freq = 326;          % Hz
stim.targetlevel = 85;    % dB SPL
stim.rampDur = 0.010;     % 10 ms ramps
stim.nTonesPerSeq = 6;
stim.nSeqPerPol = 81;
stim.interSeqSilence = 0.5; % seconds
stim.rates = [0.5 2 6];     % Hz
stim.taskDeviantProb = 0.17;
stim.deviantFreq = 380;     % Hz (task condition)



%% Create stimulus structure
stim_all = {};
stim_meta = [];

condIdx = 0;

for rate = stim.rates

    period = 1 / rate;
    stim.toneL = 0.5 * period;
    stim.silenceL = 0.5 * period;

    % generate standard & deviant tones   
    Lr = round(stim.rampDur * stim.fs);
    r = [sin(linspace(0, pi/2, Lr)), ones(1, round(stim.fs*stim.toneL) - 2*Lr), sin(linspace(pi/2, 0, Lr))];
    tone_std = r .* sin(2*pi*stim.freq * linspace(0, stim.toneL, round(stim.toneL*stim.fs)));

    tone_dev = r .* sin(2*pi*stim.deviantFreq * linspace(0, stim.toneL, round(stim.toneL*stim.fs)));

    % calibration
    choiceTT.targetlevel = stim.targetlevel + ...
        nHLfilter(stim.targetlevel, stim.freq);

    tone_std_L = calfilter('L',choiceTT.filename,choiceTT.targetlevel,tone_std,stim.fs);
    tone_std_R = calfilter('R',choiceTT.filename,choiceTT.targetlevel,tone_std,stim.fs);

    tone_dev_L = calfilter('L',choiceTT.filename,choiceTT.targetlevel,tone_dev,stim.fs);
    tone_dev_R = calfilter('R',choiceTT.filename,choiceTT.targetlevel,tone_dev,stim.fs);

    for taskFlag = [0 1]   % no-task / task
        for polarity = [1 -1]

            condIdx = condIdx + 1;
            y = [];
            trig = [];

            for seq = 1:stim.nSeqPerPol

                % decide deviant
                hasDev = taskFlag && (rand < stim.taskDeviantProb);
                if hasDev
                    devPos = randi([2 stim.nTonesPerSeq-1]);
                else
                    devPos = [];
                end

                for t = 1:stim.nTonesPerSeq

                    if hasDev && t == devPos
                        L = tone_dev_L * polarity;
                        R = tone_dev_R * polarity;
                        trigval = 200; % deviant
                    else
                        L = tone_std_L * polarity;
                        R = tone_std_R * polarity;
                        trigval = 100; % standard
                    end

                    y = [y; [L(:) R(:)]];
                    trig = [trig; trigval*ones(length(L),1)];

                    % silence after tone (except last)
                    if t < stim.nTonesPerSeq
                        ns = round(stim.silenceL * stim.fs);
                        y = [y; zeros(ns,2)];
                        trig = [trig; zeros(ns,1)];
                    end
                end

                % inter-sequence silence
                ns = round(stim.interSeqSilence * stim.fs);
                y = [y; zeros(ns,2)];
                trig = [trig; zeros(ns,1)];
            end

            stim_all{condIdx} = [y trig];
            stim_meta(condIdx).rate = rate;
            stim_meta(condIdx).task = taskFlag;
            stim_meta(condIdx).polarity = polarity;
        end
    end
end



% Plot


subplot(3,1,1)
cond1 = stim_all{:,1};
t = 0:1/stim.fs:5-1/stim.fs;
plot(t,cond1(1:5*stim.fs,1))
xlabel("Time(s)")
ylabel("Amplitude")
title("0.5 Hz stim")


subplot(3,1,2)
cond1 = stim_all{:,5};
t = 0:1/stim.fs:5-1/stim.fs;
plot(t,cond1(1:5*stim.fs,1))
xlabel("Time(s)")
ylabel("Amplitude")
title("2 Hz stim")

subplot(3,1,3)
cond1 = stim_all{:,9};
t = 0:1/stim.fs:5-1/stim.fs;
plot(t,cond1(1:5*stim.fs,1))
xlabel("Time(s)")
ylabel("Amplitude")
title("6 Hz stim")
sgtitle("5s example of sequences")


%% Playback
InitializePsychSound;
pah = PsychPortAudio('Open', [], [], 0, stim.fs, 3);

for c = 1:1 %length(stim_all)

    fprintf('Playing condition %d/%d\n',c,length(stim_all))

    PsychPortAudio('FillBuffer', pah, stim_all{c}');
    PsychPortAudio('Start', pah, 1, 0, 1);

    while PsychPortAudio('GetStatus',pah).Active
        WaitSecs(0.01);
    end

    WaitSecs(1);
end

PsychPortAudio('Close');


%% old code



% Calibration right channel
choiceTT.targetlevel = stim.targetlevel+nHLfilter(stim.targetlevel,stim.toneF); %nHL compensation (dB)
choiceTT.signal = a;
a_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
a_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);

% varying isi or 500 ms fixed?
isi = [0.5]; %isi = [.5 .6 1.4 1.5];
mean(isi)

id = 1:length(isi);
Nrep = round(stim.Lstim/(sum(isi)+length(isi)*stim.toneL));
N = Nrep*length(isi);
Lstim = (sum(isi)+length(isi)*stim.toneL)*Nrep;
Lstim = (ceil(stim.Lstim*stim.fs))/stim.fs;
fprintf('Total stimulation time: %.2f sec\n',stim.Lstim)
fprintf('Number of presentations : %d\n',N)



isi = repmat(isi,1,Nrep); % inter-stimulus-intervals
id = repmat(id,1,Nrep);
rp = randperm(length(isi));
isi = isi(rp);
id = id(rp);
ons = round([1 fs*cumsum(isi)]); % onsets in samples

% make trigger
trig_dur = round(1e-3*fs); % trigger duration in samples
trig_ampl = 1;

% make tone and trigger sequence
y = zeros((Lstim*fs),3); 
for ii = 1:length(ons)
    y(ons(ii):ons(ii)+length(a)-1,1) = a_L;
    y(ons(ii):ons(ii)+length(a)-1,2) = a_R;
    y(ons(ii):ons(ii)+trig_dur-1,3) = trig_ampl;
end
y=y';
stim.trigval = 100;
stim.isi = isi;
stim.id = id;
stim.ons = ons;
 
%% Experimental part
stim_y = y(1:2,:)';
%add trigger here
trig_y = y(3,:)';
try
    PsychPortAudio('GetOpenDeviceCount')
    PsychPortAudio('close');
catch
end
InitializePsychSound;
dev = PsychPortAudio('GetDevices');
 
%FOR PHY2
devid = 42; % 0 for PHYS2 HEAAUD
devid = 2; % my headphones mme
selectchannel =  [4 5 17;0 0 0]; %  [4 12;0 0]; ER2 + adat3
selectchannel = [0 1;0 0]; %on laptop
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
disp('Ready')
while KbCheck; end
WaitSecs(.1);
KbWait(-1);
 
Ntrials = 1;
resp_mat=[];
try
    for k=1:Ntrials
        ListenChar(2)
        % get stim from stim struct
        %stim_y = y(;
        %trig_y = stim.trigger{k};
        
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
        disp(['Playing stimulus...'])
        
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
%    stim.tones = [];
    save(['AEP_stim_',subid,'_aborted.mat'],'stim','-v7.3');
    cd ..
    psychrethrow(psychlasterror);
end
ListenChar(0);
trig.disconnect();
delete(trig);
cd('./_data')
%stim.tones = [];
save(['AEP_stim_',subid,'.mat'],'stim','-v7.3');
cd ..
disp('Experiment Done!')
%
end

