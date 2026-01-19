%% Upreat neuro pilot 1: auditory AEP sequences
% Stimulus generation and presentation script for Pilot 1: effects of rates and
% sleepiness
% plays tone sequences with rates 0.5, 2 and 6 Hz with and without deviant
% detection task


% Triggers: 
clear all
close all
clc
addpath(genpath(pwd))
 
%% General stimulus parameters (same style as yours)
stim = struct;
stim.fs = 48000;
stim.freq = 326;            % Hz (standard)
stim.targetlevel = 85;      % dB SPL
stim.rampDur = 0.010;       % 10 ms ramps
stim.nTonesPerSeq = 6;      % used for 2 & 6 Hz
stim.nSeqPerPol = 81;
stim.interSeqSilence = 0.5; % seconds (between sequences)
stim.rates = [0.5 2 6];     % Hz
stim.taskDeviantProb = 0.17;
stim.deviantFreq = 380;     % Hz (task condition)

% Behavioral
respWindow = [0, 1.2];      % seconds after deviant onset to count HIT
responseKeyName = 'SPACE';
escKeyName      = 'ESCAPE';
quitKeyName     = 'q';

% Blocks (counterbalancing)
taskFirst = false;          % true = Task first, false = No-Task first

% Calibration file
choiceTT = struct;
choiceTT.filename = 'ER3_MAG_48k_current';    % <<-- SET CAL FILE
choiceTT.targetlevel = stim.targetlevel + nHLfilter(stim.targetlevel, stim.freq);

% Audio device routing (choose your setup)
useLabDevice = true;
if useLabDevice
    devid = 42;                 % your lab device index (e.g., PHYS2)
    selectchannel = [4 5; 0 0]; % stereo only (L->4, R->5). Change if needed.
else
    devid = [];
    selectchannel = [0 1; 0 0]; % laptop default
end

%% Build stimulus structure 
stim_all = {};   % each cell: [N x 3] = [L R trig] (trig for visualization only)
stim_meta = [];  % struct per condition

condIdx = 0;

for rate = stim.rates

    period = 1 / rate;
    toneL = 0.5 * period;
    silenceL = 0.5 * period;

    % SPECIAL CASE: 0.5 Hz → single-tone "sequence"
    if abs(rate - 0.5) < 1e-6
        nTonesPerSeq_effective = 1;
    else
        nTonesPerSeq_effective = stim.nTonesPerSeq;
    end

    % generate standard & deviant tones (ramped)
    Lr = round(stim.rampDur * stim.fs);
    Ntone = round(toneL * stim.fs);
    r = [sin(linspace(0, pi/2, Lr)), ones(1, Ntone - 2*Lr), sin(linspace(pi/2, 0, Lr))];

    tvec = (0:Ntone-1)/stim.fs;
    tone_std_raw = r .* sin(2*pi*stim.freq        * tvec);
    tone_dev_raw = r .* sin(2*pi*stim.deviantFreq * tvec);

    % calibration
    tone_std_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, tone_std_raw, stim.fs);
    tone_std_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, tone_std_raw, stim.fs);

    tone_dev_L = calfilter('L', choiceTT.filename, choiceTT.targetlevel, tone_dev_raw, stim.fs);
    tone_dev_R = calfilter('R', choiceTT.filename, choiceTT.targetlevel, tone_dev_raw, stim.fs);

    for taskFlag = [0 1]   % no-task / task
        for polarity = [1 -1]

            condIdx = condIdx + 1;
            y = [];         % concatenated audio [L R]
            trig = [];      % concatenated trigger column (100/200 during tones; 0 elsewhere) — for plots
            seqOnsetsSamples = zeros(stim.nSeqPerPol,1);         % absolute seq starts in samples (within this condition)
            toneOnsetsRelPerSeq = cell(stim.nSeqPerPol,1);       % per-seq: onset samples relative to seq start
            toneCodesPerSeq     = cell(stim.nSeqPerPol,1);       % per-seq: 100/200

            cursor = 0; % absolute sample index within the condition

            for seq = 1:stim.nSeqPerPol

                % decide deviant for this sequence (task only)
                if taskFlag && (rand < stim.taskDeviantProb)
                    if nTonesPerSeq_effective == 1
                        devPos = 1;
                    else
                        devPos = randi([2 nTonesPerSeq_effective-1]); % 2..5
                    end
                    hasDev = true;
                else
                    devPos = [];
                    hasDev = false;
                end

                % store sequence onset (absolute samples)
                seqOnsetsSamples(seq) = cursor;

                toneOnsetsThis = zeros(1, nTonesPerSeq_effective);
                toneCodesThis  = zeros(1, nTonesPerSeq_effective);

                for t = 1:nTonesPerSeq_effective
                    toneOnsetsThis(t) = 0; %#ok<AGROW> % will be replaced below but keep shape
                    if hasDev && t == devPos
                        Lw = tone_dev_L(:) * polarity;
                        Rw = tone_dev_R(:) * polarity;
                        trigval = 200; % deviant
                    else
                        Lw = tone_std_L(:) * polarity;
                        Rw = tone_std_R(:) * polarity;
                        trigval = 100; % standard
                    end

                    % onset of this tone (relative to seq start)
                    if t == 1
                        relOnset = 0;
                    else
                        relOnset = (t-1) * (toneL + silenceL) - silenceL; % i.e., (t-1)*toneL + (t-2)*silenceL
                        % But with tone=silence, it's simpler: relOnset = (t-1) * (toneL + silenceL) - silenceL
                        % Actually better to accumulate:
                    end
                    % Use sample-accurate accumulation:
                    if t == 1
                        relOnsetSamples = 0;
                    else
                        relOnsetSamples = round((t-1)*(toneL + silenceL)*stim.fs) - round(silenceL*stim.fs);
                    end
                    toneOnsetsThis(t) = relOnsetSamples;

                    % append tone
                    y = [y; [Lw Rw]];
                    trig = [trig; trigval * ones(length(Lw),1)];
                    cursor = cursor + length(Lw);

                    % within-sequence silence after tone (except last)
                    if t < nTonesPerSeq_effective
                        ns = round(silenceL * stim.fs);
                        y = [y; zeros(ns,2)];
                        trig = [trig; zeros(ns,1)];
                        cursor = cursor + ns;
                    end

                    toneCodesThis(t) = trigval;
                end

                % inter-sequence silence
                nsGap = round(stim.interSeqSilence * stim.fs);
                y = [y; zeros(nsGap,2)];
                trig = [trig; zeros(nsGap,1)];
                cursor = cursor + nsGap;

                % store per-seq onsets/codes (relative to seq start)
                toneOnsetsRelPerSeq{seq} = toneOnsetsThis;
                toneCodesPerSeq{seq}     = toneCodesThis;
            end

            % store in stim_all / stim_meta (trig kept only for plotting/inspection)
            stim_all{condIdx} = [y trig];
            stim_meta(condIdx).rate      = rate;
            stim_meta(condIdx).task      = taskFlag;
            stim_meta(condIdx).polarity  = polarity;
            stim_meta(condIdx).toneL     = toneL;
            stim_meta(condIdx).silenceL  = silenceL;
            stim_meta(condIdx).nTonesPerSeq_eff = nTonesPerSeq_effective;
            stim_meta(condIdx).seqOnsetsSamples = seqOnsetsSamples;         % [nSeq x 1]
            stim_meta(condIdx).toneOnsetsRelPerSeq = toneOnsetsRelPerSeq;   % {nSeq} -> [1 x nTones]
            stim_meta(condIdx).toneCodesPerSeq     = toneCodesPerSeq;       % {nSeq} -> [1 x nTones]
        end
    end
end


%%  Plots 
% Fix: use curly braces {1}, not {:,1}
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

% Helper: Plot exactly ONE FULL sequence from each condition
plotOneSequencePerCondition(stim_all, stim_meta, stim);


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

% Assign unique SEQUENCE-ONSET TTLs per condition **within each block**
% NO-TASK: 11..16 (rate asc × polarity +1,-1)
% TASK:    21..26
seqOnsetCodes = containers.Map;  % key = sprintf('%s_%d', blockLabel, condIdx)
assignSeqOnsetCodes(blockNames, blocks, stim_meta, seqOnsetCodes);

% Initialize PsychPortAudio (stereo) and HEATriggerbox
InitializePsychSound(1);
try PsychPortAudio('Close'); catch, end
nPlayChans = size(selectchannel,2);
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
log.events   = []; % sequence onset & tone onsets
log.responses= []; % button responses (task only)
eventIdx = 0; respIdx = 0;

disp('Ready. Press any key to begin Block 1...');
while KbCheck; end
KbWait(-1);

abortAll = false;

for b = 1:2
    blockLabel = blockNames{b};
    showInstructions(blockLabel, responseKeyName);
    disp('Press ANY KEY to start this block...');
    while KbCheck; end
    KbWait(-1);

    for ii = 1:numel(blocks{b})
        c = blocks{b}(ii);
        this = stim_all{c};       % [N x 3] (L,R,trigForPlot)
        yLR  = this(:,1:2);       % play only stereo
        rate = stim_meta(c).rate;
        pol  = stim_meta(c).polarity;
        tsk  = stim_meta(c).task;

        seqOnsetsSamples  = stim_meta(c).seqOnsetsSamples;             % [nSeq x 1]
        toneOnsetsRelPerSeq = stim_meta(c).toneOnsetsRelPerSeq;        % {nSeq} -> [1 x nTones]
        toneCodesPerSeq     = stim_meta(c).toneCodesPerSeq;            % {nSeq} -> [1 x nTones]
        nSeq = numel(seqOnsetsSamples);
        nTonesPerSeq_eff = stim_meta(c).nTonesPerSeq_eff;

        % Condition-specific sequence onset TTL code
        condKey = sprintf('%s_%d', blockLabel, c);
        seqTTL  = seqOnsetCodes(condKey);

        fprintf('Block %d [%s] | Cond %d/%d: rate=%.1f, task=%d, pol=%+d | SeqOnsetTTL=%d\n', ...
            b, blockLabel, ii, numel(blocks{b}), rate, tsk, pol, seqTTL);

        % Prepare keyboard for task
        if tsk == 1
            KbQueueFlush;
            KbQueueStart;
        end

        % Start audio
        PsychPortAudio('FillBuffer', pah, yLR'); % chan x time
        [~, onsetTime] = PsychPortAudio('Start', pah, 1, 0, 1); % wait for actual start

        % Schedule TTLs: sequence-onsets + per-tone onsets
        % First tone TTL shifted by +1 ms relative to seq onset to avoid overlap.
        pulseDur = 0.001; % 1 ms

        % Precompute absolute onsets
        seqOnsetsSecAbs = onsetTime + (seqOnsetsSamples / stim.fs);
        toneOnsetsSecAbs = cell(nSeq,1);
        for s = 1:nSeq
            rel = toneOnsetsRelPerSeq{s} / stim.fs;
            % Shift first tone if its onset is == 0
            if ~isempty(rel) && rel(1) == 0
                rel(1) = rel(1) + pulseDur;
            end
            toneOnsetsSecAbs{s} = seqOnsetsSecAbs(s) + rel;
        end

        % Send sequence onset codes
        for s = 1:nSeq
            if trigbox.is_connected()
                WaitSecs('UntilTime', seqOnsetsSecAbs(s));
                trigbox.set_trigger(seqTTL);
                WaitSecs('UntilTime', seqOnsetsSecAbs(s) + pulseDur);
                trigbox.set_trigger(0);
            end
            % Log sequence-onset
            eventIdx = eventIdx + 1;
            log.events(eventIdx) = struct( ...
                'type','seqOnset','block',b,'blockLabel',blockLabel, ...
                'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                'sequence',s,'ttl',seqTTL,'time',seqOnsetsSecAbs(s));
        end

        % Send tone triggers (100 std / 200 dev)
        for s = 1:nSeq
            codes = toneCodesPerSeq{s};
            tonAbs= toneOnsetsSecAbs{s};
            for tt = 1:numel(codes)
                if trigbox.is_connected()
                    WaitSecs('UntilTime', tonAbs(tt));
                    trigbox.set_trigger(codes(tt));
                    WaitSecs('UntilTime', tonAbs(tt) + pulseDur);
                    trigbox.set_trigger(0);
                end
                % Log tone event
                eventIdx = eventIdx + 1;
                log.events(eventIdx) = struct( ...
                    'type','tone','block',b,'blockLabel',blockLabel, ...
                    'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                    'sequence',s,'toneIndex',tt,'trigCode',codes(tt),'time',tonAbs(tt));
            end
        end

        % Wait end / ESC
        sstat = PsychPortAudio('GetStatus', pah);
        while sstat.Active
            [kd,~,kc] = KbCheck(-1);
            if kd && kc(escKey)
                fprintf(2,'ESC pressed. Aborting...\n');
                abortAll = true; break;
            end
            WaitSecs(0.01);
            sstat = PsychPortAudio('GetStatus', pah);
        end
        if abortAll, break; end

        % Responses (TASK only): match to deviants in [0, 1.2] s
        if tsk == 1
            [pressed, firstPress] = KbQueueCheck;
            KbQueueStop;

            % Collect all deviant onsets for this condition
            devTimes = [];
            for s = 1:nSeq
                codes = toneCodesPerSeq{s};
                tonAbs= toneOnsetsSecAbs{s};
                dIdx = find(codes==200);
                if ~isempty(dIdx)
                    devTimes(end+1) = tonAbs(dIdx(1)); %#ok<AGROW>
                end
            end

            if pressed && firstPress(respKey) > 0
                tPress = firstPress(respKey);
            else
                tPress = [];
            end

            % Match each deviant to the first unused press in window
            used = false(size(tPress));
            for d = 1:numel(devTimes)
                if isempty(tPress), % MISS
                    respIdx = respIdx + 1;
                    log.responses(respIdx) = struct( ...
                        'block',b,'blockLabel',blockLabel,'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                        'sequence',[], 'key',responseKeyName,'tPress',NaN, ...
                        'devOnset',devTimes(d),'RT',NaN,'hit',false);
                    continue;
                end
                % Find candidate presses in [dev, dev+window]
                inWin = find(~used & tPress >= (devTimes(d)+respWindow(1)) & tPress <= (devTimes(d)+respWindow(2)), 1, 'first');
                if ~isempty(inWin)
                    RT = tPress(inWin) - devTimes(d);
                    used(inWin) = true;
                    hit = true;
                    respIdx = respIdx + 1;
                    log.responses(respIdx) = struct( ...
                        'block',b,'blockLabel',blockLabel,'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                        'sequence',[], 'key',responseKeyName,'tPress',tPress(inWin), ...
                        'devOnset',devTimes(d),'RT',RT,'hit',hit);
                else
                    % MISS
                    respIdx = respIdx + 1;
                    log.responses(respIdx) = struct( ...
                        'block',b,'blockLabel',blockLabel,'cond',c,'rate',rate,'task',tsk,'polarity',pol, ...
                        'sequence',[], 'key',responseKeyName,'tPress',NaN, ...
                        'devOnset',devTimes(d),'RT',NaN,'hit',false);
                end
            end
        end
    end

    if abortAll, break; end

    % Break between blocks
    if b == 1
        fprintf('\n=== End of Block 1 [%s] ===\n', blockLabel);
        fprintf('Take a short break. Press ANY KEY to continue or "%s" to quit.\n', quitKeyName);
        while KbCheck; end
        [~,~,kc] = KbWait(-1);
        if kc(quitKey), disp('Quit requested.'); break; end
    end
end

% Cleanup + Save
try PsychPortAudio('Stop',pah,1); catch, end
try PsychPortAudio('Close',pah); catch, end
if trigbox.is_connected(), trigbox.set_trigger(0); end
trigbox.disconnect(); delete(trigbox);
ListenChar(0);

if ~exist('./_data','dir'), mkdir('./_data'); end
save(fullfile('./_data', sprintf('AEP_%s_%s.mat', subid, datestr(now,'yyyymmdd_HHMMSS'))), ...
     'stim','choiceTT','stim_all','stim_meta','blocks','blockNames','seqOnsetCodes','log','-v7.3');

disp('Experiment Done & saved.');


%% ===== Helpers =====
function assignSeqOnsetCodes(blockNames, blocks, stim_meta, seqOnsetCodes)
% NO-TASK: 11..16 ; TASK: 21..26  (order within block: rate asc × polarity +1,-1)
for b = 1:2
    base = 10 + 10*(b-1); % 10 for NO-TASK, 20 for TASK (since blockNames defined accordingly)
    for i = 1:numel(blocks{b})
        c = blocks{b}(i);
        key = sprintf('%s_%d', blockNames{b}, c);
        seqOnsetCodes(key) = base + i; % 11..16 or 21..26
    end
end
end

function showInstructions(blockLabel, responseKeyName)
clc;
fprintf('============================================\n');
fprintf('              %s BLOCK\n', blockLabel);
fprintf('============================================\n\n');
switch upper(blockLabel)
    case 'NO-TASK'
        fprintf('- Please relax and listen passively.\n- No response is required.\n\n');
    case 'TASK'
        fprintf('- Press %s whenever you hear a DEVIANT tone.\n- Respond quickly and accurately.\n\n', responseKeyName);
end
end

function plotOneSequencePerCondition(stim_all, stim_meta, stim)
% Plot the FIRST full sequence (incl inter-seq gap) from each condition
figure('Color','w','Name','One sequence per condition');
nConds = numel(stim_all);
nRows = ceil(nConds/3);
tl = tiledlayout(nRows,3,'TileSpacing','compact','Padding','compact');

for c = 1:nConds
    rate     = stim_meta(c).rate;
    toneL    = stim_meta(c).toneL;
    silenceL = stim_meta(c).silenceL;
    nTones   = stim_meta(c).nTonesPerSeq_eff;

    % compute samples per sequence
    NsTone   = round(toneL    * stim.fs);
    NsSil    = round(silenceL * stim.fs);
    NsGap    = round(stim.interSeqSilence * stim.fs);
    NsPerSeq = nTones*NsTone + (nTones-1)*NsSil + NsGap;

    x = stim_all{c};
    if size(x,1) < NsPerSeq
        continue;
    end
    seg = x(1:NsPerSeq,:);
    t = (0:NsPerSeq-1)/stim.fs;

    nexttile;
    plot(t, seg(:,1), 'b'); hold on;
    if size(seg,2) >= 3
        yyaxis right; plot(t, seg(:,3), 'k'); ylabel('Trig');
        yyaxis left; ylabel('Amp (L)');
    else
        ylabel('Amp (L)');
    end
    title(sprintf('C%02d: rate=%.1f, task=%d, pol=%+d', c, rate, stim_meta(c).task, stim_meta(c).polarity));
    xlabel('Time (s)'); grid on;
end
sgtitle('First sequence (per condition): L audio + (optional) trig');
end
