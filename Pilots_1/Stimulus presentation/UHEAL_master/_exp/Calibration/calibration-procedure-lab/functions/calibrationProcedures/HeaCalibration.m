
classdef (Abstract) HeaCalibration < handle
    %HEACALIBRATION class (abstract)
    %
    %   HEACALIBRATION Properties:
    %       fs - Sampling frequency, in Hertz
    %       audioClass - audio utility to use (@HeaAudioPsych for example)
    %       recordIRClass - impulse response recording class (@HeaRecIRSweep for example)
    %       calculateFilterClass - inverse filter calculation class (@HeaCalcFilterMag for example)
    %
    %
    %   HEACALIBRATION Methods:
    %       calibrate - obj.calibrate() runs the whole procedure, including
    %       measuring the calibrator and saving
    %       save - obj.save() saves the object with a unique name
    %       toPDF - obj.toPDF() prints a pdf summary
    %   
    %   HEACALIBRATION Static Methods:
    %       createUserParameterFile - HeaCalibration.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %   This is the super class of several implementations, see below:
    %
    %   See also HEACALPSYBOOTH
    
    properties (Abstract)
          fs; % Sampling frequency, in Hertz
          audioClass; % audio utility to use (@HeaAudioPsych for example)
          recordIRClass; % impulse response recording class (@HeaRecIRSweep for example)
          calculateFilterClass; % inverse filter calculation class (@HeaCalcFilterMag for example)
    end
    
    properties (Constant, Access = protected)
        userParametersFileName = 'userParametersHeaCalibration.m';
        defaultParametersFileName = 'defaultParametersHeaCalibration.m';
    end
    
    methods (Abstract)
        calibrate(obj); % obj.calibrate() runs the whole procedure, including
        %measuring the calibrator and saving
        save(obj); % obj.save() saves the object with a unique name
        toPDF(obj); % obj.toPDF() prints a pdf summary
    end

    
    methods (Static, Access = protected)
        function defaultValues = loadDefaultParameters(backend)
            
            % Load default parameters
            defaultParametersHeaCalibration;
            
            % Overwrite with user parameters
            if exist('userParametersHeaCalibration','file')
                userParametersHeaCalibration;
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
            fid = fopen(HeaCalibration.defaultParametersFileName);
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
            outputFilename = HeaCalibration.userParametersFileName;
            if exist(outputFilename,'file') && ~overwrite
                error(['HeaCalibration:createUserParameterFile:The file %s ',...
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