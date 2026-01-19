% HeaCal default parameters. If you want to have different starting
% parameters, run HeaCal.createUserParameterFile and modify the file
% 'userParametersHeaCal.m' which should be now in your path.

%% General parameters

%general.fs = 48000; % Sampling frequency, in Hertz
%general.audioClass = @HeaAudioPsych; % audio utility to use (@HeaAudioPsych for example)
%general.recordIRClass = @HeaRecIRSweep; % impulse response recording class (@HeaRecIRSweep for example)
%general.calculateFilterClass = @HeaCalcFilterMag; % inverse filter calculation class (@HeaCalcFilterMag for example)

%% HeaCalPSYBooth

%PSYBooth.boothName = 'PSY1'; %ID of the booth
%PSYBooth.userEmail = 'anec@elektro.dtu.dk'; %DTU email of the person calibrating
%PSYBooth.listChannelsPlay = 1:4; %list of channels ID
%PSYBooth.listChannelsName = {{'HD650 Left' 'HD650 Right' 'HDA200 Left' 'HDA200 Right'}}; %list of channels name (warning, use two "{{")
%PSYBooth.channelRec = 1; % channels to record from
%PSYBooth.listRefSPLCalibrator = [93.8 93.8 93.8 93.8]; % list of the calibrator reference SPL values
%PSYBooth.durationCalibrator = 10; % duration of the recording of the calibrator, in seconds

%% HeaCalPHYBooth
