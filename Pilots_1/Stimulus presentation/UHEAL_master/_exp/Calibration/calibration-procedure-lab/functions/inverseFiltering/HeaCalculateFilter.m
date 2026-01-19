
classdef (Abstract) HeaCalculateFilter < handle
    %HEACALCULATEFILTER class (abstract)
    %
    %   HEACALCULATEFILTER Properties (general):
    %       fs - sampling frequency (Hz)
    %       preDelay - desired predelay before main peak in truncated impulse response 'h' (ms)
    %       decayTime - desired decay time after main peak in 'h' (ms)
    %
    %   HEACALCULATEFILTER Properties (original inpulse responses):
    %       originalImpulseResponse; % Impulse response to invert
    %       truncatedImpulseResponse; % Impulse response truncated
    %
    %   HEACALCULATEFILTER Properties (calculated):
    %       invertedImpulseResponse; % Impulse response to invert
    %       delayInverseImpulseResponse; % Delay created by the inverse response (samples)    
    %
    %   HEACALCULATEFILTER Methods:
    %       run - obj.run() runs the recording procedure for this channel
    %   
    %   HEACALCULATEFILTER Static Methods:
    %       createUserParameterFile - HeaCalculateFilter.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %   This is the super class of several implementations, see below:
    %
    %   See also HEACALCFILTERMAG
    
    properties (Abstract)
        fs; % sampling frequency (Hz)
        preDelay; % desired predelay before main peak in truncated impulse response 'h' (ms)
        decayTime; % desired decay time after main peak in 'h' (ms)
    end
    
    %Impulse responses
    properties (Abstract)
        originalImpulseResponse; % Impulse response to invert
        truncatedImpulseResponse; % Impulse response truncated
    end
    
    %Calculated filters
    properties (Abstract)
        invertedImpulseResponse; % Impulse response to invert
        delayInverseImpulseResponse; % Delay created by the inverse response (samples)
    end
    
    properties (Constant, Access = protected)
        userParametersFileName = 'userParametersHeaCalculateFilter.m';
        defaultParametersFileName = 'defaultParametersHeaCalculateFilter.m';
    end
    
    methods (Abstract)
        run(obj); % obj.run() runs the inverting procedure
    end

    
    methods (Static, Access = protected)
        function defaultValues = loadDefaultParameters(backend)
            
            % Load default parameters
            defaultParametersHeaCalculateFilter;
            
            % Overwrite with user parameters
            if exist('userParametersHeaCalculateFilter','file')
                userParametersHeaCalculateFilter;
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
            fid = fopen(HeaCalculateFilter.defaultParametersFileName);
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
            outputFilename = HeaCalculateFilter.userParametersFileName;
            if exist(outputFilename,'file') && ~overwrite
                error(['HeaCalculateFilter:createUserParameterFile:The file %s ',...
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