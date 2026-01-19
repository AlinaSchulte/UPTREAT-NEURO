
% genBMchirp.m - generates cochlea delay compensating chirp 
%
% Usage: out = out = genBMchirp(lowf, uppf, fs)
%
% lowf = lower corner frequency
% lowf = upper corner frequency
%
% out = output column vector


function out = genBMchirp(lowf, uppf, fs);

out = bmchirp(lowf,uppf,fs,0.125,1)';
