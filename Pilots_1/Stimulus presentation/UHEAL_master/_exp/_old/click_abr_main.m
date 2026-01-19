function click_abr_main(dat)
clc
 %% script to run
 %1) ABR click stim 9.1 Hz 
 %2) ABR click stim 40.1 Hz
 % 3-4 = inverted polarity

%% Subid
subid = dat.subject;
ears = {'L','R'};
choiceTT.ear = ears(dat.ear);
choiceTT.ear_noise = ears(setxor(dat.ear,[1,2]));
choiceTT.filename = 'ER3_MAG_48k_current';
 %% Generate ABR stim
 
 %stim parameters
fs = 48000;
click_dur = 0.0001; % 100 mu s
polarity = -1; %alternating polarity
rate = 9.0022; %divisable by the eeg fs %9.1; % stim presentation rate  (9.1 Hz or 40.1)
stimL = 1/rate; % stim length (s)
nreps = 3000; % repetitions (2000, 1000 in each polarity?)
peSPL = 110;%

stim = struct;

 %conditions:
 % reps x conds x Nperiods
 reps = 1;
 conds = [];
 conds = repmat([1,2,3,4],1,reps);
 %conds = Shuffle(conds);
 Nperiods = 3000;

 %%
 %tic
 % generating trials for each condition
stim = {};
stim.stimtype = [repmat({'click'},1,4)];
stim.polarity = [repmat([1,-1],1,2)];
stim.ear = repmat(dat.ear,1,4);
stim.rate = repmat([9.0022,39.9610],1,2);
stim.targetlevel = repmat(peSPL,1,4);
stim.noiselevel = repmat(55,1,4);

 for tt=1:length(conds)
     cc=conds(tt);     % condition
     trigval = conds(tt)*10+80; % Trigger value
     
     
     trig_c = [];

     %generate ABR stimuli

     if strcmp(stim.stimtype(cc),'click')
         rate = stim.rate(cc);
         stimL = 1/rate; % stim length (s)
         %gen stim
         t = 0:1/fs:stimL-1/fs;
         %click
         click =stim.polarity(cc)*[ones(round(click_dur*fs),1); ...
             zeros(round((stimL-click_dur)*fs)-1,1)];
         abr_stim = click;
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
     orig_gf = 0.5;
     %corresponded to 127.6 peSPL (dB) click
     %corresponded to 116.5 peSPL chirp
     gf_click = orig_gf*10^((peSPL-128.2)/20);
     gf_chirp = orig_gf*10^((peSPL-116.5)/20);
     if strcmp(stim.stimtype(cc),'click')
         abr_stim = abr_stim*gf_click;
     elseif strcmp(stim.stimtype(cc),'chirp')
         abr_stim = abr_stim*gf_chirp;
     end

     % Nperiods for the trial
     abr_stim = repmat(abr_stim',1,Nperiods)';
     trig_c = repmat(trig,1,Nperiods)';
     
     % add contra noise
     %TEN noise
     noiseL = length(abr_stim)/fs;
     flow = 20;fhigh=20000;
     [noise] = rand(noiseL*fs,1);%TEN_noise(noiseL,fs,flow,fhigh);
     
     % Noise Calibration chosen channel
     choiceTT.targetlevel = stim.noiselevel(cc);
     choiceTT.signal = noise;
     noise = calfilter(choiceTT.ear_noise, choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
%      % nHL compensation
%      [gainHL,f] = nHLfilter(stim.noiselevel(cc),20:12000);
%      filtercoefs = cvt_design_IG_filter(gainHL, 20:12000, fs);
%      noise = filtfilt(filtercoefs,1,noise);

     % stereo
     if stim.ear == 1
         abr_stim = [abr_stim noise];
     else
         abr_stim = [noise abr_stim];
     end
     
     % Make struct
     stim.abr_stim{tt} = abr_stim;
     stim.trigger{tt} = trig_c;
     stim.time = [0:1/fs:length(abr_stim)/fs-1/fs]';
     stim.trigval(tt) = trigval;
     stim.conds = conds;
     
 clc
 disp(num2str(tt))
     
 end
 %toc
 %reshuffling
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

%% init screen
%Screen('Preference', 'SkipSyncTests', 0);
KbName('UnifyKeyNames');
escKey = KbName('e');
quitKey = KbName('q');

%%
% ready window
disp('Ready')
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
    cd('./_data')
    stim.abr_stim = [];
    save(['click_abr_stim_',subid,'_aborted.mat'],'stim','-v7.3');
    cd ..
    PsychPortAudio('Stop', pah);
    trig.disconnect();
    delete(trig);
    while KbCheck; end
    psychrethrow(psychlasterror);
end
ListenChar(0);
trig.disconnect();
delete(trig);
cd('./_data')
stim.abr_stim = [];
save(['click_abr_stim_',subid,'.mat'],'stim','-v7.3');
cd ..
disp('Experiment Done!')
%
end
