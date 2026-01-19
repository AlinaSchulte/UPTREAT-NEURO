
classdef HeaCalcFilterLSQ < HeaCalculateFilter
    %HEACALCFILTERLSQ class
    %
    %   HEACALCFILTERLSQ Properties (general):
    %       fs - sampling frequency (Hz)
    %       preDelay - desired predelay before main peak in truncated impulse response 'h' (ms)
    %       decayTime - desired decay time after main peak in 'h' (ms)
    %
    %   HEACALCFILTERLSQ Properties (lsq specific):
    %       lowDropFreq - lower roll-off frequency (below this frequency the spectrum will be set constant)
    %       highDropFreq - upper roll-off frequency (above this frequency the spectrum will be set constant)    
    %
    %   HEACALCFILTERLSQ Properties (original inpulse responses):
    %       originalImpulseResponse; % Impulse response to invert
    %       truncatedImpulseResponse; % Impulse response truncated
    %
    %   HEACALCFILTERLSQ Properties (calculated):
    %       invertedImpulseResponse; % Impulse response to invert
    %       delayInverseImpulseResponse; % Delay created by the inverse response (samples)
    %
    %   HEACALCFILTERLSQ Methods:
    %       run - obj.run() runs the recording procedure for this channel
    %
    %   HEACALCFILTERLSQ Static Methods:
    %       createUserParameterFile - HeaCalcFilterLSQ.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %
    %   See also HEACALCULATEFILTER
    
    properties
        fs; % sampling frequency (Hz)
        preDelay; % desired predelay before main peak in truncated impulse response 'h' (ms)
        decayTime; % desired decay time after main peak in 'h' (ms)
        lowDropFreq; % lower roll-off frequency (below this frequency the spectrum will be set constant)
        highDropFreq; % upper roll-off frequency (above this frequency the spectrum will be set constant)
    end
    
    %Impulse responses
    properties
        originalImpulseResponse; % Impulse response to invert
        truncatedImpulseResponse = []; % Impulse response truncated
    end
    
    %Calculated filters
    properties
        invertedImpulseResponse = []; % Impulse response to invert
        delayInverseImpulseResponse = []; % Delay created by the inverse response (samples)
    end
    
    %Default parameters
    properties (Hidden = true)
        defaultParametersToLoad = 'lsq';
    end
    
    %HeaCalculateFilter methods
    methods
        
        function filterObj = HeaCalcFilterLSQ(varargin)
            %Constructor of the HeaRecIRSweep object
            
            % Load user/default parameters
            defaultValues = filterObj.loadDefaultParameters(filterObj.defaultParametersToLoad);
            
            %%%%%%%%%%% Check input arguments %%%%%%%%%%%%%%%%%%%%%
            p = inputParser;
            %general parameters
            addOptional(p,'fs',defaultValues.fs,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'preDelay',defaultValues.preDelay,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'decayTime',defaultValues.decayTime,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'originalImpulseResponse',defaultValues.originalImpulseResponse,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            %lsq specific
            addOptional(p,'lowDropFreq',defaultValues.lowDropFreq,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'highDropFreq',defaultValues.highDropFreq,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            
            parse(p,varargin{:})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%% Assign properties   %%%%%%%%%%%%%%%%%%%%%%
            filterObj.fs = p.Results.fs;
            filterObj.preDelay = p.Results.preDelay;
            filterObj.decayTime = p.Results.decayTime;
            filterObj.originalImpulseResponse = p.Results.originalImpulseResponse;
            filterObj.lowDropFreq = p.Results.lowDropFreq;
            filterObj.highDropFreq = p.Results.highDropFreq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function run(obj)
            %RUN obj.run() inverts
            %the impulse response
            
            % truncate the impulse response
            obj.truncatedImpulseResponse = truncir(obj.preDelay,...
                obj.decayTime,obj.originalImpulseResponse,obj.fs);
            
            % invert the impulse response
            obj.lsqinv;
            
            %limit the low and high frequencies
            obj.maglim;
            
            % calculate the group delay of the inverse filter (Re original impulse response)
            obj.delayInverseImpulseResponse = sdelay(obj.originalImpulseResponse,...
                obj.invertedImpulseResponse);
            
            % visualizing equalized responses
            visualeq(obj.originalImpulseResponse,...
                obj.invertedImpulseResponse,obj.fs);
            
        end
        
    end
    
    methods (Hidden = true)
        function lsqinv(obj)
            %LSQINV(OBJ) least squares inverse solution of an impulse response.
            %   lsqinv(h) calculates the inverse solution of the impulse response h
            %   via the pseudoinverse. The pseudoinverse provides a least squares
            %   solution.
            %
            % INPUT:
            %   truncatedImpulseResponse
            %
            % OUTPUT:
            %   invertedImpulseResponse
            %
            
            % April 2015, Henrik Gert Hassager
            % revision July 2015
            
            % create a convolution matrix for the impulse
            H = convmtx(obj.truncatedImpulseResponse(:),length(obj.truncatedImpulseResponse));
            % calculate the least squares solutions of the convolution matrix
            Hinv = (H'*H)\H';
            % calculate all solutions
            impm = H*Hinv;
            % find solution that gives the maximum impulse
            [dummy, idx] = max(diag(impm));
            % return the best inverse solution
            hinv = Hinv(:,idx);
            obj.invertedImpulseResponse = hinv(:);
        end
        
        function maglim(obj)
            %MAGLIM bandpass filters a signal with a 4th order butterworth filter and
            %   the zero-phase digital filtering function filtfilt.
            %
            % INPUT:
            %   sig     signal
            %   fcL     lower frequency for bandpass filter
            %   fcH     upper frequency for bandpass filter
            %   fs      sampling frequency
            %
            % OUTPUT:
            %   out    filtered signal
            %
            
            [b, a] =  butter(4,[obj.lowDropFreq obj.highDropFreq]/(obj.fs/2),'bandpass');
            obj.invertedImpulseResponse = filtfilt(b,a,obj.invertedImpulseResponse);
            
        end
        
    end
    
end