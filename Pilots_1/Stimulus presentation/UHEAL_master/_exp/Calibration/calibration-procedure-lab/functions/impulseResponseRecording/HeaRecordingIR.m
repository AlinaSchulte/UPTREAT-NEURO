
classdef (Abstract) HeaRecordingIR < handle
    %HEARECORDINGIR class (abstract)
    %
    %   HEARECORDINGIR Properties (general):
    %       fs - Sampling frequency, in Hertz
    %       audioClass - audio utility to use (@HeaAudioPsych for example)
    %       channelPlay - ID of the channel to calibrate
    %       channelRec - ID of the channel to record from
    %       refFreq - frequency used as a reference (Hz)
    %       refAmp - amplitude of reference signal
    %       refSPLCalibrator - calibrator output used as a reference (dB SPL)
    %       measuredRMSMic - RMS value measured by the microphone at "obj.refFreq" Hz for a dB SPL value of "obj.refSPLCalibrator"    
    %
    %   HEARECORDINGIR Properties (measured):
    %       measuredSPLChannel - dB SPL at "obj.refFreq" Hz for an rms amplitude of "obj.refAmp" of the channel being calibrated
    %       measuredImpulseResponse - Impulse response
    %
    %   HEARECORDINGIR Methods:
    %       run - obj.run() runs the recording procedure for this channel
    %   
    %   HEARECORDINGIR Static Methods:
    %       createUserParameterFile - HeaRecordingIR.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %   This is the super class of several implementations, see below:
    %
    %   See also HEARECIRSWEEP
    
    properties (Abstract)
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
    properties (Abstract)
        measuredSPLChannel; % dB SPL at "obj.refFreq" Hz for an rms amplitude of "obj.refAmp" of the channel being calibrated
        measuredImpulseResponse; % Impulse response        
    end
    
    properties (Constant, Access = protected)
        userParametersFileName = 'userParametersHeaRecordingIR.m';
        defaultParametersFileName = 'defaultParametersHeaRecordingIR.m';
    end
    
    methods (Abstract)
        run(obj); % obj.run() runs the recording procedure for this channel
    end

    
    methods (Static, Access = protected)
        function defaultValues = loadDefaultParameters(backend)
            
            % Load default parameters
            defaultParametersHeaRecordingIR;
            
            % Overwrite with user parameters
            if exist('userParametersHeaRecordingIR','file')
                userParametersHeaRecordingIR;
            end
            
            % Create the parameters structure
            generalDefaults = eval('general');
            backendDefaults = eval(backend);
            M = [fieldnames(generalDefaults)', fieldnames(backendDefaults)';
                struct2cell(generalDefaults)', struct2cell(backendDefaults)'
                ];
            defaultValues = struct(M{:});
        end
    end
    
    methods (Static)
        function createUserParameterFile(varargin)
            %CREATEUSERPARAMETERFILE  Creates a blank config file in the current folder.
            
            % Open the default parameter files
            fid = fopen(HeaRecordingIR.defaultParametersFileName);
            iLine = 1;
            tline = fgetl(fid);
            newfile = {};
            while ischar(tline)
                newfile{iLine} = regexprep(tline, '^(\w.*)', '%$1');
                tline = fgetl(fid);
                iLine = iLine + 1;
            end
            fclose(fid);
            
            % Check if the user parameters should be overwritten
            if nargin > 0
                overwrite = varargin{1};
            else
                overwrite = false;
            end
            
            
            % Write the user parameters file
            outputFilename = HeaRecordingIR.userParametersFileName;
            if exist(outputFilename,'file') && ~overwrite
                error(['HeaRecordingIR:createUserParameterFile:The file %s ',...
                    'already exists. Call the method with a "true" ', ...
                    'argument to overwrite the file.'], outputFilename);
            else
                fid = fopen(outputFilename, 'w');
            end
            
            for iLine = 1:length(newfile)
                fprintf(fid, '%s\n', newfile{iLine});
            end
            
            fclose(fid);
            
        end
        
        
        
    end
    
end