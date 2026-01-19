function [TEN_No,ratio_1k_dB] = getTEfilter(fftpts,fs,LPbin,HPbin)

% Generates filter for making threshold equalizing (TEN) noise for FTtest
% Based on TEN_cahr by Santurette/ Oxenham
% 
% IN:
% fftpts : no of fft points
% fs : sampling frequency
% LPbin : bin for lower cutoff
% HPbin : bin for high cutoff

binfactor = fftpts/fs;
f_Hz = ((1:fftpts)/binfactor)/1000; % Frequency vector in kHz

K = [0.0500   13.5000
     0.0630   10.0000
     0.0800    7.2000
     0.1000    4.9000
     0.1250    3.1000
     0.1600    1.6000
     0.2000    0.4000
     0.2500   -0.4000
     0.3150   -1.2000
     0.4000   -1.8500
     0.5000   -2.4000
     0.6300   -2.7000
     0.7500   -2.8500
     0.8000   -2.9000
     1.0000   -3.0000
     1.1000   -3.0000
     2.0000   -3.0000
     4.0000   -3.0000
     8.0000   -3.0000
    10.0000   -3.0000
    15.0000   -3.0000];

% K values are interpolated
KdB = spline(K(:,1),K(:,2),f_Hz')';

% Calculate ERB at each freq.
ERB = 24.7*((4.37*f_Hz)+1);
cr_erb = KdB+(10*log10(ERB)); % Adjustment values
TEN_No = -cr_erb; % TEN spectrum level for a 0-dB power of signal at threshold

index1kERB = f_Hz > 0.935 & f_Hz < 1.0681;  % Cams for 1 kHz +- 0.5 Cams
total_level_dB = 10.*log10(sum(10.^(TEN_No(LPbin:HPbin)/10))); % Overall level of broadband TEN noise
total_level_1k = 10.*log10(sum(10.^(TEN_No(index1kERB)/10))); % Level of TEN noise in 1-ERB band around 1-kHz
ratio_1k_dB = total_level_dB-total_level_1k; % Difference between overall level and 1-kHz band level
