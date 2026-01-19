function [irt] = truncir(t1,t2,ir,fs)
%TRUNCIR truncating a impulse response. This part is defined by the time 
%   constants t1 and t2.
%
% INPUT:
%   t1     desired time before the main peak in the IR (ms)
%   t2     desired time after the main peak in the IR (ms)
%   ir     the impulse response
%   fs     sampling frequency (Hz)
%
% OUTPUT:
%   irt    truncated impulse response
%

[m,k] = max(abs(ir));
n1 = k-round(t1*fs/1000);
if n1 < 1;
    n1 = 1;
    warning('The time constant t1 has been decreased!')
end
n2 = k+round(t2*fs/1000);
if n2 > length(ir);
    n2 = length(ir);
    warning('The time constant t2 has been decreased!')
end
irt = ir(n1:n2);