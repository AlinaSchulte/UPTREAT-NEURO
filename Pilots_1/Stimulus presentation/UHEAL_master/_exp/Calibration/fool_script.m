signal = resample(signal, 441, 480);

%% ref sig
figure(2); clf;
plot(Sig)

% load signal and look for convenient index
hold on
plot(signal)

signal = signal(1:3e5);
[xc, lag] = xcorr(signal, Sig);
[v, id] = max(abs(xc));
L = lag(id)
N = length(Sig);
subsSig = signal(L:N + L-1);

SigCpy = Sig;
Sig = subsSig;
figure(2)
hold on
plot(subsSig) 

%% sweep sig
signal = resample(signal, 441, 480);
f = figure
plot(y)
% load signal here
[xc, lag] = xcorr(signal, y);
[v, id] = max(abs(xc));
L = lag(id);
N = length(y);
subsy = signal(L:N + L-1);

yCpy = y;
y = subsy;
figure(f)
hold on
plot(subsy)

% checking sweep spectrum
a = (signal(100001:144100).*hann(44100))
freqz(a, 1, 44100, 44100)