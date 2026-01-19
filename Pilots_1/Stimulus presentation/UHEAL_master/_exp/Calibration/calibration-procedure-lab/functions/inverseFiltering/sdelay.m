function d = sdelay(h,hinv)
%SDELAY finds the delay in samples when an impulse and it's inverse are 
%   convolved.
%
% INPUT:
%   h       impulse response 
%   hinv    inverse of impulse response
%
% OUTPUT:
%   d       delay in samples
%

% delay of filter in samples
[m, d] = max(abs(conv(h,hinv)));
d = d - 1;
