function [gain,varargout] = nHLfilter(level,frequencies)

[nHL,f]=iso226(level);
%column 1 is freq
%column 2 is nHL gain

filt_c = [f;nHL-level]';

gain =spline(filt_c(:,1),filt_c(:,2),frequencies);

varargout(1) = {f};