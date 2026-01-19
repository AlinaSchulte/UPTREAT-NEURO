function AEP_EFR_main(dat)
clc
addpath(genpath(pwd))

%% script to run AEP

%8 conditions in total

%% Subid
subid = dat.subject;
%% Choose calibration file
ears = {'L','R'};
choiceTT.ear = ears(dat.ear);
choiceTT.ear_noise = ears(setxor(dat.ear,[1,2]));
choiceTT.filename = 'ER3_MAG_48k_current';

%% Generate BB stim
% stimulus parameters
stim = struct;
stim.targetlevel = 75;
stim.fs = 48000;
stim.dur = .3; % length of tones

%filter settings
[b,a] = butter(9,2000/(stim.fs/2),'low');
[bh,ah] = butter(9,500/(stim.fs/2),'high');

%% isi sweep
im = linspace(0.5,2,5);
im = flip(kron(im,ones(1,4)));
im = [im flip(im)];
im = repmat(im,1,13);
stim.isi = im;
stim.Ntrials = length(im);


%%
this_stim = [];
this_trig = [];
ons = [];

for kk = 1:stim.Ntrials
    isi_trial = stim.isi(kk);
    %% generate BB noise

    
    bbnoise = randn(stim.dur*stim.fs,1)*0.5;
    bbnoise = bbnoise';
    ramp = hann(round(0.05*stim.fs));
    win = [ramp(1:end/2);ones(length(bbnoise)-length(ramp),1);ramp(end/2+1:end)]';
    bbnoise = bbnoise.*win;
    bbnoise = bbnoise/max(abs(bbnoise))*0.5;
    
    % add modulation
    t = 0:1/stim.fs:length(bbnoise)/stim.fs-1/stim.fs;
    modl = 120;
    moddepth = 1;
    mod = (1+(moddepth)*(sin(2*pi*modl*t - pi/2)));
    mod = mod/max(abs(mod));
    bbnoise = mod.*bbnoise;
    
    bbnoise = filter(bh,ah,bbnoise);bbnoise = filter(b,a,bbnoise); % filter stimulus
    
    % Calibration
    choiceTT.targetlevel = stim.targetlevel;
    choiceTT.signal = bbnoise';
    a_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, stim.fs);
    a_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, stim.fs);
    
    bbnoise_cal = [a_L a_R]';
    
    
    % make trigger
    trig_dur = round(1e-3*stim.fs); % trigger duration in samples
    trig_stim = [ones(1,trig_dur) zeros(1,length(bbnoise)-trig_dur)];
    
    this_stim = [this_stim zeros(2,round(isi_trial*stim.fs)) bbnoise_cal];
    this_trig = [this_trig zeros(1,round(isi_trial*stim.fs)) trig_stim];
    
    ons(kk) = length([this_trig zeros(1,round(isi_trial*stim.fs))])/stim.fs; %onsets in sec
    stim.isi_cond(kk) = isi_trial;
end

%%
stim.trigval = 100;
stim.isi = im;
ids = unique(im);
id = zeros(size(stim.isi));
for cc = 1:length(ids)
    id(find(im==ids(cc))) = cc;
end
stim.id = id;
stim.ons = ons;

fprintf('Total stimulation time: %.2f min\n',length(this_stim)/stim.fs/60)
fprintf('Number of presentations : %d\n',stim.Ntrials)

%% Experimental part
stim_y = this_stim';
%add trigger here
trig_y = this_trig';
fs = stim.fs;
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
disp('Ready')
fprintf('Total stimulation time: %.2f min\n',length(this_stim)/stim.fs/60)
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

