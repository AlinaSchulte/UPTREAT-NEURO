function b = cvt_design_IG_filter(gain, aud_fc, fs)
% function b = cvt_design_IG_filter(gain, aud_fc, fs)


ntaps = 256;

aud_fc = [0 aud_fc fs/2];                           % add 0 and fs/2
aud_fc = aud_fc./(fs/2);                            % Normalize

% ODS(2, ['cbs_design: filter specs: ' num2str(AC_HL) ' dB HL']);

gain = [gain(1) gain gain(end)];                % add 0 and fs/2
gain = 10.^((gain)/20);                         % dB to power


% Design and filter
b = fir2(ntaps,aud_fc,gain);

% plot it
% if debuglevel > 2
%    figure(11);
%    freqz(b,1,512,fs);
%    set(gca, 'XScale', 'log');
%    title('IG R')
% end
