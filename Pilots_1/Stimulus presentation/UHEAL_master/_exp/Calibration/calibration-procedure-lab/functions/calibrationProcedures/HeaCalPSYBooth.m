
classdef HeaCalPSYBooth < HeaCalibration
    %HEACALPSYBOOTH class
    %
    %   HEACALPSYBOOTH Properties (General):
    %       fs - Sampling frequency, in Hertz
    %       audioClass - audio utility to use (@HeaAudioPsych for example)
    %       recordIRClass - impulse response recording class (@HeaRecIRSweep for example)
    %       calculateFilterClass - inverse filter calculation class (@HeaCalcFilterMag for example)
    %
    %   HEACALPSYBOOTH Properties (Specific to PsyBooth):
    %       boothName - ID of the booth
    %       userEmail - DTU email of the person calibrating
    %       listChannelsPlay - list of the channels to record
    %       listChannelsName - list of the name of the channels
    %       channelRec - ID of the recording channel
    %       listRefSPLCalibrator - list of the reference SPL values for the calibrator
    %       durationCalibrator - duration of the recording of the calibrator, in seconds
    %
    %   HEACALPSYBOOTH Properties (Measured):
    %       measuredRMSMic - %Measured RMS from the calibrator with the different couplers
    %       objImpResp - Impulse Response object from the obj.recordIRClass (will be the size of obj.listChannelsPlay)
    %       objInvFilter - Impulse Response object from the obj.recordIRClass (will be the size of obj.listChannelsPlay)
    %
    %   HEACALPSYBOOTH Methods:
    %       calibrate - obj.calibrate() runs the whole procedure, including
    %       measuring the calibrator and saving
    %       save - obj.save() saves the object with a unique name
    %       toPDF - obj.toPDF() prints a pdf summary
    %
    %   HEACALPSYBOOTH Static Methods:
    %       createUserParameterFile - HeaCalibration.createUserParameterFile
    %       creates a blank config file in the current folder
    %
    %   This is the super class of several implementations, see below:
    %
    %   See also HEACALIBRATION
    
    properties
        fs; % Sampling frequency, in Hertz
        audioClass; % audio utility to use (@HeaAudioPsych for example)
        recordIRClass; % impulse response recording class (@HeaRecIRSweep for example)
        calculateFilterClass; % inverse filter calculation class (@HeaCalcFilterMag for example)
        boothName; %ID of the booth
        userEmail; %DTU email of the person calibrating
        listChannelsPlay; % list of the channels to record
        listChannelsName; % list of the name of the channels
        channelRec; % ID of the recording channel
        listRefSPLCalibrator; % list of the reference SPL values for the calibrator
        durationCalibrator; % duration of the recording of the calibrator, in seconds
    end
    
    %Measured values
    properties
        measuredRMSMic = []; %Measured RMS from the calibrator with the different couplers
        objImpResp = []; % Impulse Response object from the obj.recordIRClass (will be the size of obj.listChannelsPlay)
        objInvFilter = []; % Impulse Response object from the obj.recordIRClass (will be the size of obj.listChannelsPlay)
    end
    
    %Default parameters
    properties (Hidden = true)
        defaultParametersToLoad = 'PSYBooth';
    end
    
    methods
        
        function calObj = HeaCalPSYBooth(varargin)
            %Constructor of the HeaCalPSYBooth object
            
            % Load user/default parameters
            defaultValues = calObj.loadDefaultParameters(calObj.defaultParametersToLoad);
            
            %%%%%%%%%%% Check input arguments %%%%%%%%%%%%%%%%%%%%%
            p = inputParser;
            addOptional(p,'fs',defaultValues.fs,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'audioClass',defaultValues.audioClass,@(x)validateattributes(x,...
                {'func'},{'nonempty'}))
            addOptional(p,'recordIRClass',defaultValues.recordIRClass,@(x)validateattributes(x,...
                {'func'},{'nonempty'}))
            addOptional(p,'calculateFilterClass',defaultValues.calculateFilterClass,@(x)validateattributes(x,...
                {'func'},{'nonempty'}))
            addOptional(p,'boothName',defaultValues.boothName,@(x)validateattributes(x,...
                {'char'},{'nonempty'}))
            addOptional(p,'userEmail',defaultValues.userEmail,@(x)validateattributes(x,...
                {'char'},{'nonempty'}))
            addOptional(p,'listChannelsPlay',defaultValues.listChannelsPlay,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'listChannelsName',defaultValues.listChannelsName,@(x)validateattributes(x,...
                {'char'},{'nonempty'}))
            addOptional(p,'channelRec',defaultValues.channelRec,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))
            addOptional(p,'listRefSPLCalibrator',defaultValues.listRefSPLCalibrator,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))   
            addOptional(p,'durationCalibrator',defaultValues.durationCalibrator,@(x)validateattributes(x,...
                {'numeric'},{'nonempty'}))               
            
            parse(p,varargin{:})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%% Assign properties   %%%%%%%%%%%%%%%%%%%%%%
            calObj.fs = p.Results.fs;
            calObj.audioClass = p.Results.audioClass;
            calObj.recordIRClass = p.Results.recordIRClass;
            calObj.calculateFilterClass = p.Results.calculateFilterClass;
            calObj.boothName = p.Results.boothName;
            calObj.userEmail = p.Results.userEmail;
            calObj.listChannelsPlay = p.Results.listChannelsPlay;
            calObj.listChannelsName = p.Results.listChannelsName;
            calObj.channelRec = p.Results.channelRec;
            calObj.listRefSPLCalibrator = p.Results.listRefSPLCalibrator;
            calObj.durationCalibrator = p.Results.durationCalibrator;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function calibrate(obj)
            %CALIBRATE obj.calibrate() runs the whole procedure, including
            %measuring the calibrator and saving
            
            %Measure the calibrator with the microphone for headphones (HP)
            %and earphones (EP) if needed
            
            %Set all the RMS from the calibrator to NaN
            obj.measuredRMSMic = NaN(size(obj.listChannelsPlay));
            
            %Calibrate each channel
            for idx = 1:length(obj.listChannelsPlay)
                
                %Get the reference RMS value from the calibrator
                obj.getRefRMSCalibrator(idx);
                
                % Make sure the correct channel will be measured
                hWarn = warndlg(['Measuring now: ' obj.listChannelsName{idx} ...
                    ', channel ' num2str(obj.listChannelsPlay(idx))],'!NEW CHANNEL!');
                waitfor(hWarn);
                
                %Create the impulse response object
                obj.objImpResp = [obj.objImpResp,...
                    obj.recordIRClass('fs',obj.fs,...
                    'channelPlay',obj.listChannelsPlay(idx),...
                    'channelRec',obj.channelRec,...
                    'refSPLCalibrator',obj.listRefSPLCalibrator(idx),...
                    'measuredRMSMic',obj.measuredRMSMic(idx))];
                
                %Measure the impulse response
                obj.objImpResp(idx).run;
                
                %Create the inverse filter object
                obj.objInvFilter = [obj.objInvFilter,...
                    obj.calculateFilterClass('fs',obj.fs,...
                    'originalImpulseResponse',obj.objImpResp(idx).measuredImpulseResponse)];
                
                %Calculate the inverse filter
                obj.objInvFilter(idx).run;
            end
            
        end
        
        function save(obj)
            %SAVE obj.save() saves the object with a unique name
        end
        
        function toPDF(obj)
            %TOPDF obj.toPDF() prints a pdf summary
        end
    end
    
    methods (Hidden = true)
        
        function getRefRMSCalibrator(obj,idx)
            
            % Measure calibration RMS by default
            measureRefRMSCalibrator = true;
            
            % If there is at least one measurement with the same SPL, that
            % has a reference RMS associated already: ask for using that or
            % for measuring a new one
            
            % Measurements with the same SPL profile and an existing RMS
            % value:
            measProf = find(~isnan(obj.measuredRMSMic(obj.listRefSPLCalibrator(idx) == obj.listRefSPLCalibrator)));
            if ~isempty(measProf)
                
                % Question dialog: re-measure calibration RMS?
                button = questdlg('Re-calibrate reference RMS?',...
                    'Re-measure calibration tone RMS?','Yes','No','No');
                
                % Offer options from earlier RMS values if experimenter
                % requests it
                if strcmpi(button,'no')
                    % Generate items in the drop-down list
                    listEntries = cell(1,numel(measProf));
                    for ii = measProf
                        listEntries{ii} = sprintf('%s (%3.2f RMS)',obj.listChannelsName{ii},obj.measuredRMSMic(ii));
                    end
                    
                    % List dialog box
                    [selection,ok] = listdlg('Name',sprintf('Measurements with %3.2f dB SPL reference',obj.listRefSPLCalibrator(idx)), ...
                        'ListString',listEntries,'SelectionMode','single');
                    
                    % If there was a selection: use that RMS value in
                    % current run
                    if ok
                        obj.measuredRMSMic(idx) = obj.measuredRMSMic(selection);
                        measureRefRMSCalibrator = false;
                    end
                end
            end
            
            % Measure 
            if measureRefRMSCalibrator
                obj.measuredRMSMic(idx) = measureCalibrator(obj.audioClass,...
                    obj.fs,obj.listRefSPLCalibrator(idx),obj.channelRec,obj.durationCalibrator);
            end
        end
    end
end