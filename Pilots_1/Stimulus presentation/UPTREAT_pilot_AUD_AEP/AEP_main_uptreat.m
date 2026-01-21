%% Upreat neuro pilot 1: auditory AEP sequences
% Stimulus generation and presentation script for Pilot 1: effects of rates and
% sleepiness
% plays tone sequences with rates 0.5, 2 and 6 Hz with and without deviant
% detection task


% 2 ways to send triggers:
%1) 3rd column audio channel: no distinction between conditions:
%- sequence onset: 5 ms pulse
%- standard tone 1ms pulse
% -double 1ms pulse
%2) Triggerbox: Condition with beginning of block:


%%
clear all
close all
clc
addpath(genpath(pwd))
cd('C:\GitHub\UPTREAT-NEURO\Pilots_1\Stimulus presentation\UPTREAT_pilot_AUD_AEP')

%Settings

%inLab = true;
inLab = false;

% Blocks (counterbalancing)
taskFirst = false;          % true = Task first, false = No-Task first
taskFirst = true;
%TP ID


%% General stimulus parameters (same style as yours)
stim = struct;
stim.fs = 48000;
stim.freq = 326;            % Hz
stim.targetlevel = 85;      % dB SPL
stim.rampDur = 0.010;       % 10 ms ramps
stim.nTonesPerSeq = 6;      % used for 2 & 6 Hz
stim.nSeqPerPol = 81;
stim.interSeqSilence = 0.5; % seconds (between sequences)
stim.rates = [0.5 2 6];     % Hz
stim.taskDeviantProb = 0.17;% how many deviants to be detected
stim.deviantFreq = 380;     % Hz

% Behavioral
respWindow = [0, 1.2];      % seconds after deviant onset to count HIT
responseKeyName = 'SPACE';
escKeyName      = 'ESCAPE';
quitKeyName     = 'q';


% Calibration file
choiceTT = struct;
choiceTT.filename = "ER3_MAG_48k_current.mat";
choiceTT.targetlevel = stim.targetlevel + nHLfilter(stim.targetlevel, stim.freq);

% Audio device routing
if inLab
    devid = 42;                 % lab setup
    selectchannel = [4 5; 0 0]; % present to best ear only -> 1 channel?
else
    devid = [];
    selectchannel = [0 1; 0 0]; % laptop default
end

PlayChans = size(selectchannel,2);
pah = PsychPortAudio('Open', devid, [], 0, stim.fs, nPlayChans, [], [], selectchannel);
PsychPortAudio('Volume', pah, 1.0);


%% Build stimulus structure

stim_all  = {};   % each cell: [N x 3] = [L R TriggerAudio]
stim_meta = [];   % metadata per condition
condIdx   = 0;

for rate = stim.rates

    period   = 1 / rate;
    toneL    = 0.5 * period;
    silenceL = 0.5 * period;

    % SPECIAL CASE: 0.5 Hz → single-tone "sequence"
    if abs(rate - 0.5) < 1e-6
        nTonesPerSeq_effective = 1;
    else
        nTonesPerSeq_effective = stim.nTonesPerSeq;
    end

    % Generate tones (same as before)
    Lr    = round(stim.rampDur * stim.fs);
    Ntone = round(toneL * stim.fs);
    r     = [sin(linspace(0, pi/2, Lr)), ...
        ones(1, Ntone - 2*Lr), ...
        sin(linspace(pi/2, 0, Lr))];

    tvec = (0:Ntone-1)/stim.fs;

    tone_std_raw = r .* sin(2*pi*stim.freq        * tvec);
    tone_dev_raw = r .* sin(2*pi*stim.deviantFreq * tvec);

    tone_std_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, tone_std_raw, stim.fs);
    tone_std_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, tone_std_raw, stim.fs);
    tone_dev_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, tone_dev_raw, stim.fs);
    tone_dev_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, tone_dev_raw, stim.fs);

    for taskFlag = [0 1]
        for polarity = [1 -1]

            condIdx = condIdx + 1;

            % Build audio + metadata for this condition

            y = [];  % [N x 2] audio
            seqOnsetsSamples        = zeros(stim.nSeqPerPol,1);
            toneOnsetsRelPerSeq     = cell(stim.nSeqPerPol,1);
            toneCodesPerSeq         = cell(stim.nSeqPerPol,1);

            cursor = 0;

            for seq = 1:stim.nSeqPerPol

                if taskFlag && rand < stim.taskDeviantProb
                    if nTonesPerSeq_effective == 1
                        devPos = 1;
                    else
                        devPos = randi([2 nTonesPerSeq_effective-1]);
                    end
                    hasDev = true;
                else
                    hasDev = false;
                end

                seqOnsetsSamples(seq) = cursor;
                toneOnsetsThis = zeros(1,nTonesPerSeq_effective);
                toneCodesThis  = zeros(1,nTonesPerSeq_effective);

                for t = 1:nTonesPerSeq_effective

                    if hasDev && t == devPos
                        Lw = tone_dev_L(:) * polarity;
                        Rw = tone_dev_R(:) * polarity;
                        code = 200;
                    else
                        Lw = tone_std_L(:) * polarity;
                        Rw = tone_std_R(:) * polarity;
                        code = 100;
                    end

                    if t == 1
                        relSample = 0;
                    else
                        relSample = round((t-1)*(toneL + silenceL)*stim.fs) ...
                            - round(silenceL*stim.fs);
                    end

                    toneOnsetsThis(t) = relSample;
                    toneCodesThis(t)  = code;

                    y      = [y; [Lw Rw]];
                    cursor = cursor + length(Lw);

                    if t < nTonesPerSeq_effective
                        ns = round(silenceL * stim.fs);
                        y = [y; zeros(ns,2)];
                        cursor = cursor + ns;
                    end
                end

                y = [y; zeros(round(stim.interSeqSilence*stim.fs),2)];
                cursor = size(y,1);

                toneOnsetsRelPerSeq{seq} = toneOnsetsThis;
                toneCodesPerSeq{seq}     = toneCodesThis;
            end


            % BUILD AUDIO TRIGGER CHANNEL

            fs        = stim.fs;
            Ntotal    = size(y,1);
            trigAudio = zeros(Ntotal,1);

            pulseTone = max(1, round(0.001 * fs));  % 1 ms tone pulse
            pulseSeq  = max(1, round(0.005 * fs));  % 5 ms sequence pulse
            gapPulse  = max(1, round(0.001 * fs));  % 1 ms gap (deviant double pulse)

            % -------- SEQUENCE ONSET PULSES --------
            for s = 1:numel(seqOnsetsSamples)
                idx0 = seqOnsetsSamples(s) + 1;
                idx1 = min(idx0 + pulseSeq - 1, Ntotal);
                trigAudio(idx0:idx1) = 1;
            end

            % -------- TONE ONSET PULSES --------
            for s = 1:numel(seqOnsetsSamples)

                seqAbs = seqOnsetsSamples(s);
                relOns = toneOnsetsRelPerSeq{s};
                codes  = toneCodesPerSeq{s};   % 100 = standard, 200 = deviant

                for tt = 1:numel(relOns)

                    idx0 = seqAbs + relOns(tt) + 1;
                    idx1 = min(idx0 + pulseTone - 1, Ntotal);

                    % First pulse (standard OR deviant)
                    if idx0 <= idx1
                        trigAudio(idx0:idx1) = 1;
                    end

                    % Second pulse ONLY for deviants
                    if codes(tt) == 200
                        idx0b = idx1 + gapPulse + 1;
                        idx1b = min(idx0b + pulseTone - 1, Ntotal);
                        if idx0b <= idx1b
                            trigAudio(idx0b:idx1b) = 1;
                        end
                    end

                end
            end

            % -------- STORE STIMULUS --------
            stim_all{condIdx} = [y trigAudio];


            % Store stimulus + metadata

            stim_all{condIdx} = [y trigAudio];

            stim_meta(condIdx).rate      = rate;
            stim_meta(condIdx).task      = taskFlag;
            stim_meta(condIdx).polarity  = polarity;
            stim_meta(condIdx).toneL     = toneL;
            stim_meta(condIdx).silenceL  = silenceL;
            stim_meta(condIdx).nTonesPerSeq_eff = nTonesPerSeq_effective;
            stim_meta(condIdx).seqOnsetsSamples      = seqOnsetsSamples;
            stim_meta(condIdx).toneOnsetsRelPerSeq   = toneOnsetsRelPerSeq;
            stim_meta(condIdx).toneCodesPerSeq       = toneCodesPerSeq;

        end
    end
end


%%  Plots
% stimulus examples
figure('Color','w','Name','5 s preview');
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% 0.5 Hz (cond 1): no-task,+1 (if you want exact index see order below)
nexttile;
cond1 = stim_all{1};
t = (0:size(cond1,1)-1)/stim.fs;
plot(t(1:min(end,5*stim.fs)), cond1(1:min(end,5*stim.fs),1), 'b'); grid on;
xlabel('Time (s)'); ylabel('Amp (L)'); title('0.5 Hz (first 5 s)');

% 2 Hz (~cond 5): no-task,+1 at rate=2
nexttile;
cond5 = stim_all{5};
t = (0:size(cond5,1)-1)/stim.fs;
plot(t(1:min(end,5*stim.fs)), cond5(1:min(end,5*stim.fs),1), 'r'); grid on;
xlabel('Time (s)'); ylabel('Amp (L)'); title('2 Hz (first 5 s)');

% 6 Hz (~cond 9): no-task,+1 at rate=6
nexttile;
cond9 = stim_all{9};
t = (0:size(cond9,1)-1)/stim.fs;
plot(t(1:min(end,5*stim.fs)), cond9(1:min(end,5*stim.fs),1), 'k'); grid on;
xlabel('Time (s)'); ylabel('Amp (L)'); title('6 Hz (first 5 s)');
sgtitle('5 s example of sequences');

%trigger channel
figure;
plot(trigAudio(1:min(end, round(0.5*fs))));
title('Trigger channel preview (first 0.5 s)');
xlabel('Samples'); ylabel('Trigger (0/1)');


%% Playback with blocks, HEATriggerbox TTLs, responses
subid = input('Enter subject ID (string): ','s'); if isempty(subid), subid='test'; end

% Build block orders: same condition ordering you had (rate-major: 0.5,2,6; then task=0,1; polarity=+1,-1)
% Your building loop produced this order:
%  1: 0.5Hz, no-task, +1
%  2: 0.5Hz, no-task, -1
%  3: 0.5Hz, task,    +1
%  4: 0.5Hz, task,    -1
%  5: 2Hz,   no-task, +1
%  6: 2Hz,   no-task, -1
%  7: 2Hz,   task,    +1
%  8: 2Hz,   task,    -1
%  9: 6Hz,   no-task, +1
% 10: 6Hz,   no-task, -1
% 11: 6Hz,   task,    +1
% 12: 6Hz,   task,    -1

noTaskIdx = find([stim_meta.task] == 0);
taskIdx   = find([stim_meta.task] == 1);

% sort each block by rate asc then polarity +1,-1 (for stable mapping)
[~,i1] = sortrows([[stim_meta(noTaskIdx).rate]' [stim_meta(noTaskIdx).polarity]'], [1 2]);
[~,i2] = sortrows([[stim_meta(taskIdx).rate]'   [stim_meta(taskIdx).polarity]'],   [1 2]);
noTaskIdx = noTaskIdx(i1);
taskIdx   = taskIdx(i2);

if taskFirst
    blocks = {taskIdx, noTaskIdx};
    blockNames = {'TASK','NO-TASK'};
else
    blocks = {noTaskIdx, taskIdx};
    blockNames = {'NO-TASK','TASK'};
end



% Also show the block order for clarity:
disp('=== Block orders (indices into stim_meta) ===');
for b = 1:numel(blocks)
    fprintf('%s: ', blockNames{b});
    fprintf('%d ', blocks{b});
    fprintf('\n');
end


% Initialize PsychPortAudio (stereo) and HEATriggerbox
InitializePsychSound(1);
try PsychPortAudio('Close'); catch, end

%  output 3 channels: L, R, TRIG
nPlayChans = 3;

if inLab
    devid = 42;
    selectchannel = [4 5 17; 0 0 0];
else
    devid = [];                           % default device
    selectchannel = [0 1 2; 0 0 0];       % only if system has 3 outs!
end

% Open device with 3 output channels at stim.fs
% NOTE: We don't use 'mode=1' here; leave latency setting default unless you need low‑latency.
pah = PsychPortAudio('Open', devid, [], 0, stim.fs, nPlayChans, [], [], selectchannel);
PsychPortAudio('Volume', pah, 1.0);

trigbox = HEATriggerbox();
trigbox.find_triggerbox_win();
trigbox.connect();
if ~trigbox.is_connected()
    warning('HEATriggerbox NOT connected — TTLs will not be sent.');
end

% Keyboard (KbQueue for reliable timing)
KbName('UnifyKeyNames');
respKey = KbName(responseKeyName);
escKey  = KbName(escKeyName);
quitKey = KbName(quitKeyName);
keyList = zeros(1,256); keyList(respKey)=1; keyList(escKey)=1; keyList(quitKey)=1;
KbQueueCreate(-1, keyList);

% Logs

log = struct;
log.subid = subid;
log.date  = datestr(now, 30);
log.fs    = stim.fs;
log.params= stim;
log.choiceTT = choiceTT;
log.blocks   = blockNames;
log.seqOnsetCodeMapKeys = seqOnsetCodes.keys;

% Use cell arrays to avoid type-mismatch on first assignment
log.events    = cell(0,1);   % each cell will hold one struct
log.responses = cell(0,1);

eventIdx = 0;
respIdx  = 0;


disp('Ready. Press any key to begin Block 1...');
while KbCheck; end
KbWait(-1);

abortAll = false;

for b = 1:2
    blockLabel = blockNames{b};

    % Show Instructions
    clc;
    fprintf('============================================\n');
    fprintf('              %s BLOCK\n', blockLabel);
    fprintf('============================================\n\n');
    switch upper(blockLabel)
        case 'NO-TASK'
            fprintf('- Instruct TP to relax and listen passively.\n- No response is required.\n\n');
        case 'TASK'
            fprintf('- Instruct TP to press %s whenever they hear a DEVIANT tone.\n- Respond quickly and accurately.\n\n', responseKeyName);
        otherwise
            fprintf('- Follow on-screen instructions.\n\n');
    end

    disp('Press ANY KEY to start this block...');
    while KbCheck(-1); end
    KbWait(-1);

    % --- Per-condition loop inside the block ---
    for ii = 1:numel(blocks{b})
        c = blocks{b}(ii);

        % --- Grab audio & meta
        this = stim_all{c};          % [N x 3] (L,R,trigAudio)
        yLR  = this(:,1:2);
        rate = stim_meta(c).rate;
        pol  = stim_meta(c).polarity;
        tsk  = stim_meta(c).task;

        seqOnsetsSamples    = stim_meta(c).seqOnsetsSamples;      % [nSeq x 1]
        toneOnsetsRelPerSeq = stim_meta(c).toneOnsetsRelPerSeq;   % {nSeq} -> [1 x nTones]
        toneCodesPerSeq     = stim_meta(c).toneCodesPerSeq;       % {nSeq} -> [1 x nTones]
        nSeq                = numel(seqOnsetsSamples);
        nTonesPerSeq_eff    = stim_meta(c).nTonesPerSeq_eff;

        % --- Semantic TTL (optional) ---
        if     abs(rate - 0.5) < 1e-6
            rateCode = 0;
        elseif abs(rate - 2.0) < 1e-6
            rateCode = 1;
        elseif abs(rate - 6.0) < 1e-6
            rateCode = 2;
        else
            error('Unexpected rate value: %.6f', rate);
        end
        taskCode = tsk; % 0/1
        if pol == -1
            polCode = 1;
        elseif pol == +1
            polCode = 2;
        else
            error('Unexpected polarity value: %+d', pol);
        end
        seqTTL_semantic = 100*rateCode + 10*taskCode + polCode;

        % --- Condition-specific (lookup) TTL for seq onset ---
        condKey = sprintf('%s_%d', blockLabel, c);
        seqTTL_onset = seqOnsetCodes(condKey);

        fprintf('Block %d [%s] | Cond %d/%d: rate=%.1f, task=%d, pol=%+d | SeqOnsetTTL=%d (semantic=%d)\n', ...
            b, blockLabel, ii, numel(blocks{b}), rate, tsk, pol, seqTTL_onset, seqTTL_semantic);

        % --- Prepare keyboard for task ---
        if tsk == 1
            KbQueueFlush;
            KbQueueStart;
        end

        % --- Start audio (use device channel count you opened with) ---
        buf  = this(:, 1:nPlayChans);
        PsychPortAudio('FillBuffer', pah, buf');   % chan x time
        [onsetTime] = PsychPortAudio('Start', pah, 1, 0, 1);

        % --- Precompute absolute seq/tone onsets (aligned to onsetTime) ---
        pulseDur = 0.001; % 1 ms
        fs = stim.fs;     % Ensure stim.fs is defined
        seqOnsetsSecAbs = onsetTime + (seqOnsetsSamples / fs);
        toneOnsetsSecAbs = cell(nSeq,1);
        for s = 1:nSeq
            rel = toneOnsetsRelPerSeq{s} / fs;
            if ~isempty(rel) && rel(1) == 0
                rel(1) = rel(1) + pulseDur;  % avoid overlap with seq TTL
            end
            toneOnsetsSecAbs{s} = seqOnsetsSecAbs(s) + rel;
        end

        % --- Send sequence-onset TTLs and log (log INSIDE the loop to keep s/t_fire in scope) ---
        if trigbox.is_connected()
            for s = 1:nSeq
                t_fire = seqOnsetsSecAbs(s);
                WaitSecs('UntilTime', t_fire);
                trigbox.set_trigger(seqTTL_onset);
                WaitSecs('UntilTime', t_fire + pulseDur);
                trigbox.set_trigger(0);

                % Log per-sequence onset
                eventIdx = eventIdx + 1;
                log.events{eventIdx} = struct( ...
                    'type','seqOnset','block',b,'blockLabel',blockLabel, ...
                    'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                    'sequence',s,'ttl',seqTTL_onset,'time',t_fire);
            end
        else
            % If no trigbox: still log expected times (optional)
            for s = 1:nSeq
                eventIdx = eventIdx + 1;
                log.events{eventIdx} = struct( ...
                    'type','seqOnset','block',b,'blockLabel',blockLabel, ...
                    'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                    'sequence',s,'ttl',seqTTL_onset,'time',seqOnsetsSecAbs(s));
            end
        end

        % --- Watch for ESC/q during playback ---
        abortedThisCond = false;
        while true
            sstat = PsychPortAudio('GetStatus', pah);
            if ~sstat.Active
                break; % playback ended normally
            end

            [keyIsDown, ~, keyCode] = KbCheck(-1);
            if keyIsDown
                if keyCode(KbName('ESCAPE')) || keyCode(KbName('q'))
                    PsychPortAudio('Stop', pah, 1);   % stop now and wait
                    abortedThisCond = true;
                    abortAll = true;
                    break;
                end
                % Drain key state
                while KbCheck(-1); end
            end
            WaitSecs(0.005);
        end
        if abortedThisCond
            break;  % leave the ii-loop; outer loop will see abortAll
        end

        % --- Ensure playback is fully stopped ---
        sstat = PsychPortAudio('GetStatus', pah);
        while sstat.Active
            [kd,~,kc] = KbCheck(-1);
            if kd && kc(KbName('ESCAPE'))
                fprintf(2,'ESC pressed. Aborting...\n');
                abortAll = true; break;
            end
            WaitSecs(0.01);
            sstat = PsychPortAudio('GetStatus', pah);
        end
        if abortAll
            break; % leave ii-loop
        end

        % --- Responses (TASK only): match to deviants in window ---
        if tsk == 1
            [pressed, firstPress] = KbQueueCheck;
            KbQueueStop;

            % Collect first deviant onset per seq (code==200)
            devTimes = [];
            for s = 1:nSeq
                codes = toneCodesPerSeq{s};
                tonAbs= toneOnsetsSecAbs{s};
                dIdx = find(codes==200, 1, 'first');
                if ~isempty(dIdx)
                    devTimes(end+1) = tonAbs(dIdx); %#ok<AGROW>
                end
            end

            % Extract press times for the response key
            if pressed && exist('respKey','var') && firstPress(respKey) > 0
                tPress = firstPress(respKey);
            else
                tPress = [];
            end

            % Match each deviant to the first unused press in [dev + respWindow(1), dev + respWindow(2)]
            used = false(size(tPress));
            for d = 1:numel(devTimes)
                if isempty(tPress)
                    % MISS
                    respIdx = respIdx + 1;
                    log.responses{respIdx} = struct( ...
                        'block',b,'blockLabel',blockLabel,'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                        'sequence',[], 'key',responseKeyName,'tPress',NaN, ...
                        'devOnset',devTimes(d),'RT',NaN,'hit',false);
                    continue;
                end
                inWin = find(~used & tPress >= (devTimes(d)+respWindow(1)) & tPress <= (devTimes(d)+respWindow(2)), 1, 'first');
                if ~isempty(inWin)
                    RT = tPress(inWin) - devTimes(d);
                    used(inWin) = true;
                    respIdx = respIdx + 1;
                    log.responses{respIdx} = struct( ...
                        'block',b,'blockLabel',blockLabel,'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                        'sequence',[], 'key',responseKeyName,'tPress',tPress(inWin), ...
                        'devOnset',devTimes(d),'RT',RT,'hit',true);
                else
                    % MISS
                    respIdx = respIdx + 1;
                    log.responses{respIdx} = struct( ...
                        'block',b,'blockLabel',blockLabel,'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                        'sequence',[], 'key',responseKeyName,'tPress',NaN, ...
                        'devOnset',devTimes(d),'RT',NaN,'hit',false);
                end
            end
        end

    end % ---------- end for ii ----------

    if abortAll
        break; % leave block loop
    end

    % --- Break between blocks: only after Block 1 ---
    if b == 1
        fprintf('\n=== End of Block 1 [%s] ===\n', blockLabel);
        fprintf('Take a short break. Press ANY KEY to continue or "%s" to quit.\n', quitKeyName);

        % Clear any residual key state
        while KbCheck(-1); end

        % Wait for any key, then check which
        KbWait(-1);
        [~,~,kc] = KbCheck(-1);

        if kc(quitKey)
            disp('Quit requested.');
            abortAll = true;
            break;
        end
        if kc(KbName('ESCAPE'))
            fprintf(2,'ESC pressed at block break. Aborting...\n');
            abortAll = true;
            break;
        end
    end

end % ---------- end for b ----------


%% Cleanup + Save
try PsychPortAudio('Stop',pah,1); catch, end
try PsychPortAudio('Close',pah); catch, end
if trigbox.is_connected(), trigbox.set_trigger(0); end
trigbox.disconnect(); delete(trigbox);
ListenChar(0);

if ~exist('./_data','dir'), mkdir('./_data'); end
save(fullfile('./_data', sprintf('AEP_%s_%s.mat', subid, datestr(now,'yyyymmdd_HHMMSS'))), ...
    'stim','choiceTT','stim_all','stim_meta','blocks','blockNames','seqOnsetCodes','log','-v7.3');

disp('Experiment Done & saved.');


