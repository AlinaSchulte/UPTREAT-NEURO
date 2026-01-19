
classdef HeaRecIRSweep < HeaRecordingIR
    %HEARECIRSWEEP HeaRecordingIR class based on sweep recording
    %
    %   HEARECIRSWEEP Properties (general):
    %       fs - Sampling frequency, in Hertz
    %       audioClass - audio utility to use (@HeaAudioPsych for example)
    %       channelPlay - ID of the channel to calibrate
    %       channelRec - ID of the channel to record from
    %       refFreq - frequency used as a reference (Hz)
    %       refAmp - amplitude of reference signal
    %       refSPLCalibrator - calibrator output used as a reference (dB SPL)
    %       measuredRMSMic - RMS value measured by the microphone at "obj.refFreq" Hz for a dB SPL value of "obj.refSPLCalibrator"    
    %
    %   HEARECIRSWEEP Properties (measured):
    %       measuredSPLChannel - dB SPL at "obj.refFreq" Hz for an rms value of 1 of the channel being calibrated
    %       measuredImpulseResponse - Impulse response
    %
    %   HEARECIRSWEEP Properties (specific to the sweep method):
    %       durationRecording -  recording duration for headphone output (s)
    %       nRepSweep - Number of sweep repetitions
    %       nPlace - number of placements for each headphone channel
    %       sweepLen - length of log sweep (s)
    %       scale - normalize impulse response? {'yes' | 'no'}
    %       preDelay - desired predelay before main peak in truncated impulse response 'h' (ms)
    %       decay - desired decay time after main peak in 'h' (ms)
    %
    %   HEARECIRSWEEP Methods:
    %       run - obj.run() runs the recording procedure for this channel
    %
    %   HEARECIRSWEEP Static Methods:
    %       createUserParameterFile - HeaRecIRSweep.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %
    %   See also 
    
    properties
        fs; % Sampling frequency, in Hertz
        audioClass; % audio utility to use (@HeaAudioPsych for example)
        channelPlay; % ID of the channel to calibrate
        channelRec; % ID of the channel to record from
        refFreq; % frequency used as a reference (Hz)
        refAmp; % amplitude of reference signal
        refSPLCalibrator; % calibrator output used as a reference (dB SPL)
        measuredRMSMic; % RMS value measured by the microphone at "obj.refFreq" Hz for a dB SPL value of "obj.refSPLCalibrator"
    end
    
    %measured properties
    properties
        measuredSPLChannel = []; % dB SPL at "obj.refFreq" Hz for an rms value of 1 of the channel being calibrated
        measuredImpulseResponse = []; % Impulse response
    end
    
    %sweep specific properties
    properties
        durationRecording;  % recording duration for headphone output (s)
        nRepSweep; % Number of sweep repetitions
        nPlace; % number of placements for each headphone channel
        sweepLen; % length of log sweep (s)
        scale; % normalize impulse response? {'yes' | 'no'}
        preDelay; % desired predelay before main peak in truncated impulse response 'h' (ms)
        decay; % desired decay time after main peak in 'h' (ms)
    end
    
    %Default parameters
    properties (Hidden = true)
        defaultParametersToLoad = 'sweep';
    end
    
    %HeaCal methods
    methods
       
        function calObj = HeaRecIRSweep(varargin)
            %Constructor of the HeaRecIRSweep object
            
            % Load user/default parameters
            defaultValues = calObj.loadDefaultParameters(calObj.defaultParametersToLoad);
            
            %%%%%%%%%%% Check input arguments %%%%%%%%%%%%%%%%%%%%%
            p = inputParser;
            %general parameters
            addOptional(p,'fs',defaultValues.fs,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'channelPlay',defaultValues.channelPlay,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'channelRec',defaultValues.channelRec,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'audioClass',defaultValues.audioClass,@(x)validateattributes(x,...
                {'func'},{'nonempty'}))
            addOptional(p,'refFreq',defaultValues.refFreq,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'refAmp',defaultValues.refAmp,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'refSPLCalibrator',defaultValues.refSPLCalibrator,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'measuredRMSMic',defaultValues.measuredRMSMic)
            %sweep specific
            addOptional(p,'durationRecording',defaultValues.durationRecording,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'nRepSweep',defaultValues.nRepSweep,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'nPlace',defaultValues.nPlace,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'sweepLen',defaultValues.sweepLen,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'scale',defaultValues.scale,@(x)validateattributes(x,...
                {'char'},{'nonempty'}))
            addOptional(p,'preDelay',defaultValues.preDelay,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'decay',defaultValues.decay,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            
            parse(p,varargin{:})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%% Assign properties   %%%%%%%%%%%%%%%%%%%%%%
            calObj.fs = p.Results.fs;
            calObj.channelPlay = p.Results.channelPlay;
            calObj.channelRec = p.Results.channelRec;
            calObj.audioClass = p.Results.audioClass;
            calObj.refFreq = p.Results.refFreq
            calObj.refAmp = p.Results.refAmp
            calObj.refSPLCalibrator = p.Results.refSPLCalibrator;
            calObj.measuredRMSMic = p.Results.measuredRMSMic;
            calObj.durationRecording = p.Results.durationRecording;            
            calObj.nRepSweep = p.Results.nRepSweep;
            calObj.nPlace = p.Results.nPlace;
            calObj.sweepLen = p.Results.sweepLen;
            calObj.scale = p.Results.scale;
            calObj.preDelay = p.Results.preDelay;
            calObj.decay = p.Results.decay;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function run(obj)
            %RUN obj.run() measures
            %the impulse response
            
            % Create reference signal
            refSig = obj.refAmp*sin((0:1:round(obj.durationRecording*obj.fs)-1)/obj.fs*2*pi*obj.refFreq)';
            % Measure headphone frequency response
            HPSig = zeros((obj.durationRecording-1)*obj.fs,obj.nPlace);
            rmsHPSig = zeros(1,obj.nPlace);
            
            H = zeros(ceil((obj.preDelay + obj.decay)/1000*obj.fs/2),obj.nPlace);
            h = zeros(ceil((obj.preDelay + obj.decay)/1000*obj.fs),obj.nPlace);
            hWarn = warndlg('Place headphone on coupler.','PLACE HEADPHONE');
            waitfor(hWarn);
            
            % Create the audio player
            player = obj.audioClass('channelsPlay',obj.channelPlay,'channelsRec',obj.channelRec,'fs',obj.fs);
            
            for ll = 1:obj.nPlace
                
                clc; fprintf('Measuring Channel %s \n',obj.channelPlay);
                
                % Play and record reference signal for  "obj.nPlace" times
                Sig = player.playrec(refSig);
                HPSig(:,ll) = Sig(obj.fs/2+1:(obj.durationRecording-1)*obj.fs+obj.fs/2);
                
                % Calculate root mean square value of recorded signal
                rmsHPSig(:,ll) = rms(HPSig(:,ll));
                
                % Measure impulse response
                [IR(:,ll),x,y,t,obj.nRepSweep,D,obj.fs] = measIR(obj,obj.nRepSweep,obj.sweepLen,obj.fs,obj.scale);
               
                %Produces a one-sided amplitude spectrum from an impulse response
                [H(:,ll),h(:,ll),fNorm] = obj.getH(obj.preDelay,obj.decay,IR(:,ll),obj.fs,obj.refFreq);
                
                if ll < obj.nPlace
                    clc; fprintf('Measurement %d/%d done, please reposition earphone, \nthen press any key to start next measurement\n',ll,obj.nPlace);
                    pause;
                elseif ll == obj.nPlace
                    clc; fprintf(['SPL of ' num2str(obj.refFreq) ' Hz sinusoid with RMS of 1:\n']);
                end
                
            end
            
            % Calculate output level of RMS 1 sinusoid (measuredSPLChannel)
            rmsHPSigmean = mean(rmsHPSig);
            obj.measuredSPLChannel = obj.refSPLCalibrator+20*log10(sqrt(2)*rmsHPSigmean/(obj.refAmp*obj.measuredRMSMic));
            disp(obj.measuredSPLChannel)
            
            % Calculate mean values of TFs for calibrating Channel
            HCh = mean(abs(H(:,1:obj.nPlace)),2);      % mean (measured) Magnitude of left EP
            obj.measuredImpulseResponse = mean(h(:,1:obj.nPlace),2);  % (truncated) impulse response of left EP
             
        end
        
    end
    
    methods (Hidden = true, Access = protected)
       
        function [IR,x,y,t,nrep,D,fs] = measIR(obj,nrep,D,fs,scale)
            % Measurement of impulse responses.
            %
            % Input
            % nrep          repetition of sweeps
            % D             approx. duration of IR in seconds (D should be larger)
            % fs            sampling frequency (Hz)
            % scale         scaling of impulse responses {'yes' | 'no'}
            
            % Output
            % IR         measured impulse response
            % x          input signal (sweep presented via soundcard output)
            % y          measured signal (measurement of the sweep played back via headphone)
            % t          time vector (for plotting)
            % nrep       number of sweep presentations
            % D          duration of each sweep (seconds)
            % fs         sampling frequency
            %
            
            
            gdB = -30;  %attenuation of the sound signal in dB (gdB > 0 clipping effects occur)
            
            nfft = D*fs;
            out = 1;
            while out < nfft;
                out = out*2;
            end
            nfft = out;
            
            D0 = .1; % silence interval (seconds) before and after the input signal
            z0 = zeros(round(D0*fs),1);
            [x, InvRef] = obj.makesweep(nfft, nrep, fs);
            x = [z0 ; 10.^(gdB/20)*x ; z0];
            N = length(x);
            
            % Create the audio player
            player = obj.audioClass('channelsPlay',obj.channelPlay,'channelsRec',obj.channelRec,'fs',obj.fs);
            
            % Play and record reference signal 
            y = player.playrec(x);
            
            M = y(:,1);
            M(N) = 0;
            IR = obj.sweep2imp(M, InvRef, nrep);
            
            if strcmpi(scale,'yes')
                IR = .99*IR/max(max(abs(IR)));
            end
            t = linspace(0,D,D*fs); % time vector [s]
            
        end
        
        function [s, H] = makesweep(obj, nfft, nrep, fs)
            % deleted input parameters: fstart, fstop, timew, amplw, end_delay, verbose
            %Creates a logarithmic (or modified logarithmic) sweep excitation signal for
            %impulse measurement purposes. Frequency domain sweep generation is used to allow
            %more control on the sweep properties. As a default the time domain sweep will
            %have a (nearly) constant amplitude.
            %
            %Using a logarithmic sweep allows the separation of distortion components - created
            %for example by a loudspeaker - from the linear response of the device under test
            %(DUT). If the sweep is long enough compared to the response of the DUT, the
            %distortion components will manifest themselves as separate components in the end
            %of the impulse response recovered by using SWEEP2IMP.
            %
            %The function allows modifications of the sweep rate as a function of frequency.
            %With such modifications the amount of energy fed to the DUT, and thus effectively
            %the signal-to-noise ratio (SNR), can be controlled as a function of frequency.
            %However, the user should be aware that such processing will change the behavior
            %of the distortion components. A single harmonic distortion component will no
            %longer start from a single instant. Instead the start time will depend on the
            %sweep rate at each frequency.
            %
            %Return values:
            %  s = The sweep sequence.
            %
            %  H = Inverse reference spectrum needed for recovering the impulse response
            %    from the measurement results. See SWEEP2IMP.
            %
            %Basic parameters:
            %  nfft = Number of samples in one sweep. For faster operation use a power of two.
            %    Nfft must be higher than the length of the impulse response being measured
            %    in order to avoid circular wrapping of the response when averaging results of
            %    several sweep repetitions.
            %
            %  nrep = Number of sweep repetitions. When averaging data from several repetitions
            %    using SWEEP2IMP, the first sequence is by default discarded.
            %
            %  fs = Sampling frequency in Hz. Default: 48000.
            %
            %Band limited sweeps:
            %  The parameters fstart and fstop allow limiting the frequency range of the sweep
            %  excitation. The sweeps always start at near DC and extend to near Nyquist
            %  frequency in order to avoid abrupt starting and stop phenomena. However, the
            %  sweeping up until fstart and from fstop on can be made fast so that little
            %  time (signal energy) is consumed at uninteresting frequencies.
            %
            %  These band limiting features will manifest themselves as zero phase band pass
            %  filtering on the final recovered impulse response. Both highpass and lowpass
            %  filters have a frequency response corresponding to a second order butterworth
            %  filter. For very narrow bands and low upper limits the band limiting may cause
            %  artifacts in the beginning and in the end of the recovered impulse response
            %  due to periodic frequency domain processing.
            %
            %  fstart = -3 dB cutoff frequency of the optional highpass filter. Default:
            %    10 Hz (corresponds to about -0.25 dB @ 20 Hz).
            %
            %  fstop = -3 dB cutoff frequency of the optional lowpass filter. Default:
            %    fs/2 (off).
            %
            %Arbitrary sweep rate and amplitude control:
            %  The parameters timew and amplw allow frequency dependent weighting of the
            %  excitation signal magnitude spectrum. These weighting functions should be
            %  column vectors of nfft+1 elements such that the first element corresponds to
            %  DC and the last to Nyquist frequency (sorry for the weird parameter form).
            %  Both weightings will be compensated when recovering the impulse response
            %  with SWEEP2IMP. See also SWEEPWEIGHT.
            %
            %  timew = Sweep magnitude spectrum weighting that will affect the sweep rate
            %    of a constant amplitude sweep.
            %
            %  amplw = Sweep magnitude spectrum weighting applied as amplitude control of
            %    the sweep signal.
            %
            %The rest:
            %  end_delay = Defines the fraction of zero samples in the end of a single
            %    period of a sweep signal. The zeros are needed to capture the full response
            %    when a single sweep is used. When averaging repeated sweeps it is best to
            %    keep this parameter zero. Defaults: when nrep = 1: 0.4, otherwise: 0.
            %
            %  verbose = Controls plotting of additional figures for debugging purposes.
            %    Value 1 can be used to investigate the effects of band limiting and user
            %    defined spectral weighting.
            %
            %References: Müller & Massarani: "Transfer-Function Measurements with Sweeps,"
            %  J. Audio Eng. Soc., Vol. 49, No. 6, Pp. 443-471, June 2001.
            %
            %See also: SWEEP2IMP, SWEEPWEIGHT
            %
            %Juha Merimaa 03.06.2003
            
            ignore_bp_in_ref = 1;
            
            % deleted input parameters are defined here:
            fstart = 10;
            fstop = fs / 2;
            timew = 0;
            amplw = 0;
            verbose = 0;
            
            if nrep > 1
                end_delay = 0;
            else
                end_delay = 0.4;
            end
            
            % allow some time for fluctuations before and after the actual sweep within the
            % final excitation signal
            fstart = max(fstart, fs / (2 * nfft));
            predelay = ceil(max(fs / fstart, nfft / 200));
            if predelay > nfft / 10
                predelay = ceil(nfft/10);
            end
            postdelay = ceil(max(fs / fstop, nfft / 200));
            if postdelay > nfft / 10
                postdelay = ceil(nfft/10);
            end
            % the number of samples in the actual sweep
            nsweep = round((1 - end_delay) * nfft) - predelay - postdelay;
            
            % transform into seconds
            sweep_sec = nsweep / fs;
            predelay_sec = predelay / fs;
            
            % from here on use a double length time window to prevent time domain artifacts
            % from wrapping to the sweep signal note: this makes nfft always even!
            nhalffft = nfft + 1; % up to and including Nyquist
            nfft = 2 * nfft;
            f = (0:nhalffft-1)' * fs / nfft;
            
            % construct a pink magnitude spectrum
            H = [fstart / f(2); sqrt(fstart ./ f(2:nhalffft))];
            
            % frequency control options that will shape the time evolution of the sweep:
            % highpass filtering
            if fstart > fs / nfft
                [bhp, ahp] = butter(2, 2*fstart/fs, 'high');
                H = H .* abs(freqz(bhp, ahp, 2*pi*f/fs));
            end
            % lowpass filtering
            if fstop < fs / 2
                [blp, alp] = butter(2, 2*fstop/fs);
                H = H .* abs(freqz(blp, alp, 2*pi*f/fs));
            end
            % user defined spectral weighting
            if timew
                H = H .* timew(:);
            end
            
            % calculate group delay
            C = sweep_sec ./ sum(H.^2);
            tg = C * cumsum(H.^2);
            tg = tg + predelay_sec;
            if verbose
                vfig = figure;
                subplot(2, 2, 1), semilogx(f(1:nhalffft), tg);
                title('Constructed group delay');
                xlabel('Frequency / Hz');
                ylabel('Time / s');
                axis tight;
                grid on;
            end
            % calculate phase
            ph = -2*pi*fs/nfft * cumsum(tg);
            % force the phase to zero at Nyquist
            ph = ph - f/f(nhalffft) .* mod(ph(nhalffft), 2*pi);
            clear tg;
            
            % optional spectral weighting controlling the amplitude of the sweep
            if amplw
                H = H .* amplw(:);
            end
            
            % create double-sided spectrum
            H = H .* exp(j*ph);
            H(nhalffft+1:nfft) = conj(H(nhalffft-1:-1:2));
            
            % convert to time domain
            s = real(ifft(H));
            
            if verbose > 1
                figure;
                subplot(2, 1, 1), plot(s, 'b-')
                title('The sweep signal before and after time windowing');
                xlabel('Time / samples');
                ylabel('Amplitude');
                hold on;
                subplot(2, 1, 2), semilogx(f(2:nhalffft), 20*log10(abs(H(2:nhalffft))+realmin), 'b-');
                xlabel('Frequency / Hz');
                ylabel('Magnitude / dB');
                axis tight;
                grid on;
                hold on;
            end
            
            % window the fluctuations before and after the actual sweep
            w = hann(2 * predelay);
            s(1:predelay) = s(1:predelay) .* w(1:predelay);
            stopind = nsweep + predelay;
            w = hann(2 * postdelay);
            s(stopind+1:stopind+postdelay) = s(stopind+1:stopind+postdelay)...
                .* w(postdelay+1:2*postdelay);
            s(stopind+postdelay+1:nfft) = 0;
            
            if verbose > 1
                subplot(2, 1, 1), plot(s, 'r-')
                hold off;
                temp = fft(s);
                temp = temp(1:nhalffft);
                subplot(2, 1, 2), semilogx(f, 20*log10(abs(temp)), 'r-');
                clear temp
                hold off;
            end
            
            % back to the user defined nfft
            nfft = nfft / 2;
            s = s(1:nfft);
            
            % normalize the amplitude and create the repetions
            normfact = 1.02 * max(abs(s));
            s = s ./ normfact;
            if nrep > 1
                s2 = s;
                for k=2:nrep
                    s2 = [s2; s];
                end
                s = s2;
                clear s2;
            end
            
            H = 1 ./ fft(s(1:nfft));
            f = (0:nfft-1)' * fs / nfft;
            % calculate a new reference spectrum without the bandpass filtering
            % (in order not to amplify noise by the inverse of the sweep energy far
            % outside the sweep band)
            if ignore_bp_in_ref
                if fstart > fs / nfft
                    H = H .* abs(freqz(bhp, ahp, 2*pi*f/fs));
                end
                if fstop < fs / 2
                    H = H .* abs(freqz(blp, alp, 2*pi*f/fs));
                end
            end
            
            if verbose
                figure(vfig);
                subplot(2, 2, 2), plot(s(1:nfft));
                title('Sweep signal, single repetition');
                xlabel('Time / samples');
                ylabel('Amplitude');
                
                temp = fft(s(1:nfft)) .* H;
                subplot(2, 2, 4), semilogx(f, 20*log10(abs(temp)+realmin));
                title('Magnitude spectrum of the recovered impulse response');
                xlabel('Frequency / Hz');
                ylabel('Magnitude / dB');
                axis tight;
                a = axis;
                a(2) = fs / 2;
                a(3) = max(a(3), -60);
                axis(a);
                grid on;
                
                temp = real(ifft(temp));
                subplot(2, 2, 3), plot(-20:20, temp([nfft-19:nfft 1:21]));
                title('Recovered impulse response');
                xlabel('Time / samples');
                ylabel('Amplitude');
                
                
            end
        end
        
        function impresp = sweep2imp(obj, s, InvRef, nrep)
            
            %Recover an impulse response from a sweep excitation measurement.
            
            % Input:
            % s        Response to the sweep signal.
            % InvRef   Inverse reference spectrum, given by MAKESWEEP.
            % nrep     Number of repetitions in the sweep excitation.
            
            %Juha Merimaa 03.06.2003
            
            if nrep > 1;
                nfft = length(InvRef);
                s2 = 1 / (nrep - 1) * s(nfft+1:2*nfft);
                for k=2:nrep-1;
                    s2 = s2 + 1 / (nrep-1) * s(k*nfft+1:(k+1)*nfft);
                end
                s = s2;
            end
            
            impresp = real(ifft(fft(s) .* InvRef));
        end
        
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
    end
    
    
end