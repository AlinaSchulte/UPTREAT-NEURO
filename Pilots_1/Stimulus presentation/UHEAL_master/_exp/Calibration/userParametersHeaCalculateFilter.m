% HeaCal default parameters. If you want to have different starting
% parameters, run HeaRecordingIR.createUserParameterFile and modify the file
% 'userParametersHeaRecordingIR.m' which should be now in your path.

%% General parameters

general.fs = 20000; % Sampling frequency, in Hertz
general.preDelay = 2; % desired predelay before main peak in truncated impulse response 'h' (ms)
general.decayTime = 98; % desired decay time after main peak in 'h' (ms)
general.originalImpulseResponse = []; % Impulse response

%% Magnitude inversion parameters

mag.filterLen = 2048; % desired length of inverse filter (samples)
mag.lowDropFreq = 50; % lower roll-off frequency (below this frequency the spectrum will be set constant)
mag.highDropFreq = 9998; % upper roll-off frequency (above this frequency the spectrum will be set constant)

%% LSQ (least squares, time domain)

lsq.lowDropFreq = 30; % lower roll-off frequency (below this frequency the spectrum will be set constant)
lsq.highDropFreq = 9998; % lower roll-off frequency (below this frequency the spectrum will be set constant)
