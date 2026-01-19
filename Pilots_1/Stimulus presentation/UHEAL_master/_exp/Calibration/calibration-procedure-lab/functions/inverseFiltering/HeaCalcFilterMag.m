
classdef HeaCalcFilterMag < HeaCalculateFilter
    %HEACALCFILTERMAG class
    %
    %   HEACALCFILTERMAG Properties (general):
    %       fs - sampling frequency (Hz)
    %       preDelay - desired predelay before main peak in truncated impulse response 'h' (ms)
    %       decayTime - desired decay time after main peak in 'h' (ms)
    %     
    %   HEACALCFILTERMAG Properties (mag specific):
    %       filterLen - desired length of inverse filter (samples)
    %       lowDropFreq - lower roll-off frequency (below this frequency the spectrum will be set constant)
    %       highDropFreq - upper roll-off frequency (above this frequency the spectrum will be set constant)
    %
    %   HEACALCFILTERMAG Properties (original inpulse responses):
    %       originalImpulseResponse; % Impulse response to invert
    %       truncatedImpulseResponse; % Impulse response truncated
    %
    %   HEACALCFILTERMAG Properties (calculated):
    %       invertedImpulseResponse; % Impulse response to invert
    %       delayInverseImpulseResponse; % Delay created by the inverse response (samples)    
    %
    %   HEACALCFILTERMAG Methods:
    %       run - obj.run() runs the recording procedure for this channel
    %   
    %   HEACALCFILTERMAG Static Methods:
    %       createUserParameterFile - HeaCalcFilterMag.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %
    %   See also HEACALCULATEFILTER
    
    properties
        fs; % sampling frequency (Hz)
        preDelay; % desired predelay before main peak in truncated impulse response 'h' (ms)
        decayTime; % desired decay time after main peak in 'h' (ms)
        filterLen; % desired length of inverse filter (samples)
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
        defaultParametersToLoad = 'mag';
    end
    
    %HeaCalculateFilter methods
    methods
       
        function filterObj = HeaCalcFilterMag(varargin)
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
            %mag specific
            addOptional(p,'filterLen',defaultValues.filterLen,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'lowDropFreq',defaultValues.lowDropFreq,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'highDropFreq',defaultValues.highDropFreq,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))            
            addOptional(p,'originalImpulseResponse',defaultValues.originalImpulseResponse,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))               
            
            parse(p,varargin{:})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%% Assign properties   %%%%%%%%%%%%%%%%%%%%%%
            filterObj.fs = p.Results.fs;
            filterObj.preDelay = p.Results.preDelay;
            filterObj.decayTime = p.Results.decayTime;
            filterObj.filterLen = p.Results.filterLen;
            filterObj.lowDropFreq = p.Results.lowDropFreq;
            filterObj.highDropFreq = p.Results.highDropFreq;
            filterObj.originalImpulseResponse = p.Results.originalImpulseResponse;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function run(obj)
            %RUN obj.run() inverts
            %the impulse response
            
            % truncate the impulse response
            obj.truncatedImpulseResponse = truncir(obj.preDelay,...
                obj.decayTime,obj.originalImpulseResponse,obj.fs);
            
            % invert the impulse response
            obj.maginv;
            
            % calculate the group delay of the inverse filter (Re original impulse response)
            obj.delayInverseImpulseResponse = sdelay(obj.originalImpulseResponse,...
                obj.invertedImpulseResponse);
            
            % visualizing equalized responses
            visualeq(obj.originalImpulseResponse,...
                obj.invertedImpulseResponse,obj.fs);
            
        end
        
    end
    
    methods (Hidden = true)
        function maginv(obj)
            %MAGINV(OBJ) inverse magnitude solution of an impulse response.
            %   maginv calculates the inverse magnitude solution of the impulse
            %   response h in the frequency domain.
            %
            % INPUT:
            %   truncatedImpulseResponse
            %   lowDropFreq     lower roll-off frequency (below this frequency the
            %                   spectrum will be set constant)
            %   highDropFreq    upper roll-off frequency (above this frequency the
            %                   spectrum will be set constant)
            %   fs              sampling frequency (Hz)
            %   filterLen       desired length of inverse filter (samples)
            %
            % OUTPUT:
            %   invertedImpulseResponse
            %
            
            n = length(obj.truncatedImpulseResponse);
            H = fft(obj.truncatedImpulseResponse,n);
            % Number of unique frequency bins
            K = (n+bitand(n,1))/2+~bitand(n,1);
            % Return positive frequencies
            H = abs(H(1:K));
            % remove the frequency roll-off for very low and very high frequencies
            N = size(H,1);
            n1 = round(obj.lowDropFreq/obj.fs*2*N+1);
            n2 = round(obj.highDropFreq/obj.fs*2*N+1);
            H(1:n1) = H(n1);
            H(n2:N) = H(n2);
            % normalized frequency
            fNorm = linspace(0,1,length(H))';
            % calculate inverse filter
            Bi = fir2(obj.filterLen,fNorm,1./H);
            % minimum phase
            [r,hinv] = rceps(Bi);
            hinv = hinv(1:obj.filterLen);
            obj.invertedImpulseResponse = hinv(:);
        end
    end
    
end