        function [H,h,fNorm] = getH(obj, t1,t2,IR,fs,refFreq)
            
            %Produces a one-sided amplitude spectrum from an impulse
            %response. A part of the impulse response is only used in the spectrum
            %derivation. This part is defined by the time constants t1 and t2.
            %
            % Input:
            % t1     Desired time before the main peak in the IR (ms)
            % t2     Desired time after the main peak in the IR (ms)
            % IR     the (measured) impulse response
            % fs     sampling frequency (Hz)
            %
            % Output
            % H      (half-sided) FFT spectrum (from 0Hz to fs/2)
            % h      truncated (and actually employed) impulse response
            % fNorm      frequency vector of amplitude spectrum 0..1 (where 1 refers to
            %       fs/2)
            
            % Update {jc} : H is now the one-sided COMPLEX VALUED FFT, not only absolute value.
            % Correction of roll-off and smoothing now moved to inverse filtering process in
            % 'HP_design_inv_Filt', 'HP_deconv_wav'
            
            % Update {gl} : H is scaled to yield 0 dB amplification at refFreq
            
            %Cutting the "important" part out of the impulse response
            
            [m,k] = max(abs(IR));
            n1 = k-round(t1*fs/1000);
            if n1 < 1;
                n1 = 1;
                warning('The time constant t1 has been decreased!')
            end
            n2 = k+round(t2*fs/1000);
            if n2 > length(IR);
                n2 = length(IR);
                warning('The time constant t2 has been decreased!')
            end
            h = IR(n1:n2);
            
            %Producing an impulse response with an even number of samples
            if rem(length(h),2) ~= 0;
                h = h(1:end-1);
            end
            
            N = length(h)/2;
            
            % frequency vector from 0...1, as needed for filter design later on
            fNorm = linspace(0,1,N)';
            
            %getting one-sided FFT spectrum
            H = fft(h);
            % LAST MODIFICATION #1: THESE LINES ADDED
            v = abs(fNorm*fs/2-refFreq);
            H = H / abs(H(find(v==min(v)))); % normalizing H
            h = ifft(H);
            % MODIFICITION #1 END
            H = H(1:N);
            
        end