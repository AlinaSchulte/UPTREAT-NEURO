% HeaCal default parameters. If you want to have different starting
% parameters, run HeaRecordingIR.createUserParameterFile and modify the file
% 'userParametersHeaRecordingIR.m' which should be now in your path.

global choice   % required to set reference frequencz from GUI (see lines 18 ff)

%% General parameters

general.fs = 48000; % Sampling frequency, in Hertz
general.channelPlay = 1; % ID of the channel to calibrate
general.channelRec = 1; % ID of the channel to record from
general.audioClass = @HeaAudioPsych; % audio utility to use (from the 'HeaAudio' class)
general.refFreq = choice.refFreq; % used to be 1000% frequency used as a reference (Hz)
general.refAmp = 0.1; % amplitude of reference signal
general.refSPLCalibrator = 93.8; % calibrator output used as a reference (dB SPL)
general.measuredRMSMic = []; % dB SPL at "obj.refFreq Hz" for an rms amplitude of "obj.refRMS" of the microphone

%% Additional lines to incorporate variable calibrator frequency (anec)
% (Onli if GUI has been called previously)
if exist('choice','var'),
    
    if isfield(choice,'refFreq');
    
        general.refFreq     = choice.refFreq;
        
    end
    
end


%% HeaVocSweep

sweep.durationRecording = 3; % recording duration for headphone output (s)
sweep.nRepSweep = 5; % Number of sweep repetitions
sweep.nPlace = 3; % number of placements for each headphone channel
sweep.sweepLen = 1; % length of log sweep (s)
sweep.scale = 'yes'; % normalize impulse response? {'yes' | 'no'}
sweep.preDelay = 2; % desired predelay before main peak in truncated impulse response 'h' (ms)
sweep.decay = 98; % desired decay time after main peak in 'h' (ms)
