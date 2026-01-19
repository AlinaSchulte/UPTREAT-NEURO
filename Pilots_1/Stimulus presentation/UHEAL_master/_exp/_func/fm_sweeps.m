function [x,sweepf] = fm_sweeps(dur,fs,flow,fhigh,dir)
%% fm-sweep stim
%fm_sweeps(dur,fs,flow,fhigh,dir)
%dur = length of sweep
%dir = direction (1)=low to high (2)=high to low

% from Krishnana&Parkinson2000
% dur 80ms with 10ms ramp in/out
% sweep 400-600Hz
% cosine-squared gating

t = 0:1/fs:dur-1/fs;

if dir
sweepf = linspace(flow,fhigh,length(t));
else
sweepf = linspace(fhigh,flow,length(t));   
end

modl=1/dur;
tonel = length(t);
r = nramp(modl,fs,0.025,tonel/fs)';
r(1) = 0;
%cosgate = (cos(2*pi*1/dur/2*t).^120-1)*-1;
gatedur = 5e-3;
gatet = 0:1/fs:gatedur*2;
cosgate = abs(cos(2*pi*1/gatedur/4*gatet).^2-1);
coswin =[cosgate(1:end/2) ones(length(t)-2*length(cosgate(1:end/2)),1)' cosgate(end/2+1:end)];

stim = coswin.*cos(2*pi.*cumsum(sweepf)/fs);
x = repmat(stim,1,8);

%plot(stim)
%hold on
%plot(coswin)
end
%sound(stim,fs)



