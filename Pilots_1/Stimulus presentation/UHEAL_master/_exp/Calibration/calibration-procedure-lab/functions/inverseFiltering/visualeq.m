function visualeq(ir,irinv,fs)

ch = conv(ir,irinv);

L = numel(ch);
NFFT = 2^nextpow2(L);
H = fft(ch,NFFT);
f = fs/2*linspace(0,1,NFFT/2+1);
Hpos = H(1:NFFT/2+1);

figure
subplot(311)
mag = 20*log10(abs(Hpos));
semilogx(f,mag)
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
ylim([-5 5])
xlim([20 fs/2])
grid on
% remove linear phase delay   
[m,k] = max(abs(ch));
w = 2*pi*f(1:NFFT/2+1)*(k-1);
ph = angle(Hpos.*exp(-1i*w/fs)');
subplot(312)
semilogx(f,ph)
xlabel('Frequency [Hz]')
ylabel('Phase [radian]')
ylimits = get(gca,'Ylim');
ydiff = ylimits(2) - ylimits(1);
if ydiff < pi + 0.2;
    ylimits(2)  = ylimits(2) + pi + 0.1;
    ylimits(1)  = ylimits(1) - pi - 0.1;
end
ylim(ylimits)
xlim([20 fs/2])
grid on
gd = -gradient(unwrap(angle(Hpos)))./(2*pi*gradient(f)')*1000;
subplot(313)
semilogx(f,gd)
xlabel('Frequency [Hz]')
ylabel('Group Delay [msec]')
ylimits = get(gca,'Ylim');
ydiff = ylimits(2) - ylimits(1);
if ydiff < 10;
    ylimits(2)  = ylimits(2) + 5;
    ylimits(1)  = ylimits(1) - 5;
end
if ydiff > 100;
    ylimits(2)  = median(gd) + 50;
    ylimits(1)  = median(gd) - 50;
end
ylim(ylimits)
xlim([20 fs/2])
grid on