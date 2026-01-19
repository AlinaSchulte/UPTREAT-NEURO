function AEP_main(dat)
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
 
%% Generate tone stim
% stimulus parameters
stim = struct;
stim.targetlevel = 70;
stim.fs = 48000;
stim.Lstim = 240; % length of stimulation
stim.toneL = .1; % length of tones
stim.toneF = 1000; % tone freq
Lr = round(.01*stim.fs); % length of tone ramp
fprintf('Total stimulation time: %d sec\n',stim.Lstim)

fs = 48000;
% generate one tone
r = [sin(linspace(0, pi/2, Lr)) ones(1,stim.fs*stim.toneL-Lr*2) ,...
    sin(linspace(pi/2, 0, Lr))];


a = r.*sin(2*pi*stim.toneF*linspace(0,stim.toneL,stim.toneL*stim.fs));


% Calibration right channel
choiceTT.targetlevel = stim.targetlevel+nHLfilter(stim.targetlevel,stim.toneF); %nHL compensation (dB)
choiceTT.signal = a;
a_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
a_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);


isi = [.5 .6 1.4 1.5];
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

