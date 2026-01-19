function click_abr_demo(dat)
clc

%This function is called by exp_dorun.m to create the 9.1 and 40.1 stimulus
%, condition order, and manage the presentation of the click abr stimulus.

%UPDATE: SAM 16.12.2019
%Added true alternating polarity ('polarity') variable was repurposed to
%toggle this feature on.
%Added 'slow attenuator' mimic feature of the eclipse whereby the stimulus
%fades in over 10 seconds before triggers begin.
%Added optional stimulus jitter feature which injects delay into the
%stimulus and associated trigger rate.

%UPDATE: SAM 19/06/2020
%Added slow ramp attemnuator for end of click train as well as the begining


%% Subid
subid = dat.subject;
%% Choose calibration file
ears = {'L','R'};
choiceTT.ear = ears(dat.ear);
choiceTT.ear_noise = ears(dat.ear);
choiceTT.filename = 'ER3_MAG_48k_current';
%% Controls 
%stim parameters
fs = 48000; 
click_dur = 0.0001; % 100 mu s
polarity = -1; %1= pos only, -1 = alt polarity
jitter = 1; % 1 = on

nreps = 500; % repetitions total
peSPL = 115.5;%

stim = struct;

%conditions:
% reps x conds x Nperiods
reps = 1; %number of repititions
conds = [];
conds = repmat([1,2],1,reps);
%conds = Shuffle(conds); %randomise order


%% generating trials for each condition
stim = {};
stim.stimtype = [repmat({'click'},1,3)];

if polarity == -1
stim.polarity = [repmat([-1,-1 -1],1,reps)]; %indicates true alt pol
elseif polarity == 1
stim.polarity = [repmat([1,1,1],1,reps)]; %only positive
else
    fprintf('Polarity variable badly entered')
    return
end

stim.ear = repmat(dat.ear,1,3*reps);
stim.rate = repmat([9.1022,40.0587],1,reps);
stim.targetlevel = repmat(peSPL,1,3*reps);
stim.noise = repmat([0 0],1,reps);
stim.noiselevel = 70;

 for tt=1:length(conds)
     cc=conds(tt);     % condition
     trigval = conds(tt)*10+40; % Trigger value + 40 (from ffr)    
     trig_c = [];

     %generate ABR stimuli

     if strcmp(stim.stimtype(cc),'click')
         rate = stim.rate(cc);
         stimL = 1/rate; % stim length (s)
         %gen stim
         t = 0:1/fs:stimL-1/fs;
         
         % generate one click epoch and trigger
         click =stim.polarity(cc)*[ones(round(click_dur*fs),1); ...
             zeros(round((stimL-click_dur)*fs)-1,1)];
         abr_stim = -click;
         trig =[ones(1,round(1e-3 * fs)) zeros(1,length(abr_stim)-round(1e-3 * fs))];
         
     elseif strcmp(stim.stimtype(cc),'chirp')        
         rate =stim.rate(cc);
         stimL = 1/rate; % stim length (s)
         %gen stim
         chp = genBMchirp(450,10000,fs);
         chp = chp(1:round(stimL*fs)-1);
         chp = stim.polarity(cc)*[chp;zeros(round(stimL-length(chp)/fs)*fs,1)];
         abr_stim = chp;
         trig =[ones(1,round(1e-3 * fs)) zeros(1,length(abr_stim)-round(1e-3 * fs))];
     else
         error('no matching stim type')
     end

     % gain correction (calibration)
     if strcmp(choiceTT.ear,'R') % right ear
     orig_gf =0.34; %check click_calib 
     else
     orig_gf = 0.286;   % left ear
     end
     
     if strcmp(stim.stimtype(cc),'click')
         abr_stim = abr_stim*orig_gf;%gf_click;
     elseif strcmp(stim.stimtype(cc),'chirp')
         abr_stim = abr_stim*gf_chirp;
     end

     %save a copy of a single epoch
     epoch = abr_stim; 
     
     %% generate whole stimulus to full length
     
     jits = zeros(nreps-1,1); %default to no jitters
     if jitter == 1
     rng('default') %reset random number generation seed for replicable results
     jits = randi(144,[nreps-1,1]); %uniform distribution of jitters 1 to 144 samples (3ms)
     end
     
    trig_c = trig; %initiate a copy of trigger pulse 
    if polarity == -1 %true alternating
        for n = 1:(nreps-1)
           epoch = -epoch;
           abr_stim = [abr_stim;zeros(jits(n),1);epoch];
           trig_c = [trig_c,zeros(1,jits(n)),trig];
        end 
       
    elseif polarity == 1 %positive only
        for n = 1:(nreps-1)
           abr_stim = [abr_stim;zeros(jits(n),1);epoch];
           trig_c = [trig_c,zeros(1,jits(n)),trig];
        end
    end
    %% generate noise
    
    noiseL = ceil(length(abr_stim)/fs);
    flow = 20;fhigh=20000;
    [noise,~,~] = TEN_noise(noiseL,fs,flow,fhigh);
    choiceTT.signal = noise;
    noise = calfilter(choiceTT.ear, choiceTT.filename, stim.noiselevel, choiceTT.signal, fs);

    %% Add slow attenuator
    ramp_reps = ceil(10*fs/(2*length(epoch))); %repetitions of two epochs (+, -)
    ramp_reps = 2*(round(ramp_reps/2)); %make sure it's an even number
    ramp = linspace(0,1,ramp_reps*2*length(epoch))'.*repmat([-epoch;epoch],ramp_reps,1);
    ramp_noise = linspace(0,1,ramp_reps*2*length(epoch))'.*noise(1:length(ramp));
    
    %prepend it to the stimulus
    abr_stim = [ramp;abr_stim];
    trig_c = [zeros(1,length(ramp)),trig_c]; %extend back trig channel with zeros
    noise = [ramp_noise;noise];
    
    %append it to the stimulus (flipped to be a fade OUT)
    abr_stim = [abr_stim; flip(ramp)];
    trig_c = [trig_c,zeros(1,length(ramp))]; %extend forwards trig channel with zeros
    noise = [noise;flip(ramp_noise)];
    
     %% stereo
     if stim.ear == 1
         abr_stim = [abr_stim zeros(size(abr_stim))];
     elseif stim.ear == 2
         abr_stim = [zeros(size(abr_stim)) abr_stim];
     end
     %% add noise
     if stim.noise(cc)
         if stim.ear(cc) == 1
             noise = [noise zeros(size(noise))];
             abr_stim = abr_stim+noise(1:length(abr_stim),:);
         elseif stim.ear(cc) ==2
             noise = [zeros(size(noise)) noise];
             abr_stim = abr_stim+noise(1:length(abr_stim),:);            
         end
     end
     
     %% Make struct
     stim.abr_stim{tt} = abr_stim;
     stim.trigger{tt} = trig_c';
     stim.time = [0:1/fs:length(abr_stim)/fs-1/fs]';
     stim.trigval(tt) = trigval;
     stim.conds = conds;
     
 clc
 disp(num2str(tt))
     
 end

%% reshuffling
%stim.stimtype = stim.stimtype{}; % OBS if multiple stimtypes
stim.polarity = stim.polarity(conds);
stim.ear = stim.ear(conds);
stim.targetlevel = stim.targetlevel(conds);
stim.trigval = stim.trigval';
stim.dat = dat;
stim.experiment_time = dat.experiment_time;


 %%
 ltmp = [];
 for ii = 1:length(stim.abr_stim)
     ltmp(ii) = length(stim.abr_stim{ii})/fs;
%      soundsc(stim.ABRstim{ii},fs)
%      pause()
 end
 clc
 disp(['length of experiment: ',num2str(sum(ltmp)/60),' min'])

 
 
 %% Experimental part
stim_y = stim.abr_stim{1};
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
selectchannel =  [4 5 0;0 0 0]; %  [4 12;0 0]; ER2 + adat3
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

%% init screen
%Screen('Preference', 'SkipSyncTests', 0);
KbName('UnifyKeyNames');
escKey = KbName('e');
quitKey = KbName('q');

%%
% ready window
clc
disp(['length of experiment: ',num2str(sum(ltmp)/60),' min'])

while KbCheck; end
WaitSecs(.1);
KbWait(-1);

Ntrials = length(stim.abr_stim);
try
    for k=1:Ntrials
        ListenChar(2)
        % get stim from stim struct
        stim_y = stim.abr_stim{k};
        trig_y = stim.trigger{k};
        
         % set trigger value
        if trig.is_connected()
          trig.set_trigger(stim.trigval(k));
        end
%         
        %fill the buffer
        PsychPortAudio('FillBuffer', pah,[stim_y trig_y]')
        %PsychPortAudio('FillBuffer', pah,[stim_y'])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Drawing screen
        clc
        disp(['Playing stimulus condition ' stim.stimtype{k} ', ' num2str(stim.rate(k)) ' Hz rate' ', Polarity ' num2str(stim.polarity(k))])

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
    ListenChar(0);
    %cd('./_data')
    stim.abr_stim = [];
    %dont save demo
    %save(['click_abr_stim_',subid,'_aborted.mat'],'stim','-v7.3');
    %cd ..
    PsychPortAudio('Stop', pah);
    %trig.disconnect();
    %delete(trig);
    while KbCheck; end
    psychrethrow(psychlasterror);
end
ListenChar(0);
%trig.disconnect();
%delete(trig);
%cd('./_data')
stim.abr_stim = [];
% dont save demo
%save(['click_abr_stim_',subid,'.mat'],'stim','-v7.3');
%cd ..
disp('Sound demo done!')
%
end
