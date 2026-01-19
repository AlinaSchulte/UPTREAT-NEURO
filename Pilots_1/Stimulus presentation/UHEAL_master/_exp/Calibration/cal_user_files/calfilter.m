function [output_signal, output_channel] = calfilter(ear, headphones_name, target_level, input_signal, fs)

% calfilter.m - filters and scales a signal to equalize the headphone
% response and set the SPL to a specifid value
% INPUT:
% ear               - the headphone side. Possible values: right: 'r', 'R', left: 'l', 'L'
% headphones_name   - name of the file containing calibration data, e.g. HD650_MAG_48k_current
% target_level      - target level fo the equalized signal (in SPL)
% input_signal      - needs to e a single column-vector
% fs                - sampling freq of the input signal (needs to match filter's fs)
% OUTPUT:
% output_signal
% output_channel    
% Author:
% Michal Fereczkowski, 24 Apr 2017, based on earlier versions from CAHR
% Edited:
% Andreas Eckey, 29 June 2017,  adding output channel specification
%                               changing string reading for loading calObj
%                               

%% check input_signal
% signal needs to be stored in a single column
if ~iscolumn(input_signal)
    if size(input_signal,1)==1
        input_signal = input_signal';
%        warning('The input signal needs to be a single column - transpositioning')
    else
        error('The input signal needs to be a single column')
    end
end

%% set RMS to 1
input_signal = input_signal/rms(input_signal);

%% load headphone equalization file and check for matching fs
% exemplary headphones_name = 'HD650_MAG_48k_current.mat';
load([headphones_name]);
if fs ~= calObjMag_48k.fs
    error('sampling frequencies of the signal ad the filter do not match!')
end

%% retrieve calibration levels and equalization filters
callevelL = calObjMag_48k.objImpResp(1).measuredSPLChannel; % Headphone calibration level at 1 kHz (dB SPL) (LEFT)
callevelR = calObjMag_48k.objImpResp(2).measuredSPLChannel; % Headphone calibration level at 1 kHz (dB SPL) (RIGHT)

eqL = calObjMag_48k.objInvFilter(1).invertedImpulseResponse; % Headphone equalization filter (LEFT)
eqR = calObjMag_48k.objInvFilter(2).invertedImpulseResponse; % Headphone equalization filter (RIGHT)

% choose side
if strcmpi(ear, 'l')
  eqF               = eqL;
  callevel          = callevelL;
  output_channel    = calObjMag_48k.listChannelsPlay(1);
elseif strcmpi(ear, 'r')
  eqF               = eqR;
  callevel          = callevelR;
  output_channel    = calObjMag_48k.listChannelsPlay(2);
else
  error('Unexpected ear string')
end

%% filter the signal with the equalization filter
filtered_signal = filter(eqF,1,input_signal);

%% set the signal level to the target value
output_signal = filtered_signal * 10^((target_level - callevel)/20);

