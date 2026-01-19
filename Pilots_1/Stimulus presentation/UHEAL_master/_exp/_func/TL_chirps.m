
function stim = TL_chirps(dur,fs,fl,fu,len)


%rising chirp from Troels Lindgren
%ltfatstart
%%
close all
%fs = 44100;
%fu = 1400;
%fl = 300;
%dur = .5;
t = 0:1/fs:len-1/fs;


phi1 = (((fu+fl).*t)./2);
phi2= ((dur.*(fu-fl).*sin((2*pi.*t)/dur))./8*pi)*4/3;


gatedur = 5e-3;
gatet = 0:1/fs:gatedur*2;
cosgate = abs(cos(2*pi*1/gatedur/4*gatet).^2-1);
coswin =[cosgate(1:end/2) ones(length(t)-2*length(cosgate(1:end/2)),1)' cosgate(end/2+1:end)];

stim = coswin.*sin(2*pi.*phi1-phi2);%cumsum(phi));
%figure(1)
%plot(t,stim)
%plot(phi1)
%hold on
%plot(phi2)
%plot(phi1-phi2+3/4*pi)


%plot(t,phi/fs)
%figure(2)
%sound(stim*0.5,fs)
%sgram(stim,fs,'wlen',round(15e-3*fs),'dynrange',50)
%ylim([100 1600])
end