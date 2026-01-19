function out = nramp(f,fs,lr,varargin)

% Creates cosine on-off ramps periodically repeating at frequency f
%
% INPUTS:
%   f   - repetition frequency in Hz
%   fs  - sampling frequency
%   lr  - length of ramp in sec
%   L   - length of output in sec
%           if default/unspecified then the output is some whole no of periods
%           and can be repeated without sampling errors
%
% jensh-0206-2019

if lr>(.5/f)
    error('length of ramp less than half a period')
end

nn = find(rem([1:500]*(fs/f),1)==0,1,'first'); % n whole periods in fs_o
n = find(rem([1:1e4]./f,1)==0,1,'first'); % some fs where we can make one period with whole n samples

for ii=1:10000
    if n>10000 % n should be at least this fs
        break
    end
    n = n*ii;
end

pl = n/f; % the length of one period
lr = round(lr*n); % length of ramp in samples
r = [sin(linspace(0, pi/2, lr)) ones(1,pl-lr*2) sin(linspace(pi/2, 0, lr))]'; % one period
r = repmat(r,nn,1); % repeat one period to match n whole periods in fs
out = resample(r,fs,n);
out = out-min(out); out = out/max(out); % scale to 0-1

if nargin>3
    L = varargin{1};
    N = ceil(L*fs);
    reps = ceil(N/length(out));
    out = repmat(out,reps,1);
    out = out(1:N);
end