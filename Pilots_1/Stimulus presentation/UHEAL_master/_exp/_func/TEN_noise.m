function [noise,noise_rms,ratio_1k] = TEN_noise(noiseL,fs,flow,fhigh)


noiseL = noiseL*fs; % noise length


lowcut = round(flow*(noiseL/fs))+1;  %low bin to cut from
highcut = round(fhigh*(noiseL/fs))+1;
bw_bins = (highcut-lowcut)+1; % no of bins in pass band

Y = zeros(1,noiseL); % create noise in spectral domain
Y(:,lowcut:highcut) = randn(1,bw_bins) + 1i*randn(1,bw_bins); % random amplitudes and phases

[TEfilter,ratio_1k] = getTEfilter(noiseL,fs,lowcut,highcut);
Y(lowcut:highcut) = Y(lowcut:highcut).*10.^(TEfilter(lowcut:highcut)/20); % threshold equalization of spectrum

noise = real(ifft(Y)); % time domain noise signal
noise_rms = sqrt(sum(noise .^ 2) ./ noiseL); % rms without gap

%soundsc(noise,fs)
end