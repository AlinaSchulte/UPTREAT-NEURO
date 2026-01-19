%% Preference
% Calibrating
addpath(genpath(pwd))
% Name to save          

fileIDstring            = sprintf('%s%scalibration_tracker%s%s_%s_MAG_48k_%s.mat',choice.path,filesep,filesep,choice.booth,choice.phone,datestr(now,'yymmdd_HHMM'));
fileIDstring_current    = sprintf('%s%scal_user_files%s%s_MAG_48k_current.mat',choice.path,filesep,filesep,choice.phone);
if choice.fs ~= 48000
    fs_txt = num2str(choice.fs);
    fileIDstring  = strrep(fileIDstring, '48000', fs_txt);
    fileIDstring_current  = strrep(fileIDstring_current, '48000', fs_txt);
end

filename                = fileIDstring;
filename_current        = fileIDstring_current;

% Create object (using defaultParametersHeaCalibration.m)
calObjMag_48k           = HeaCalPSYBooth;

%% Object will now be modified to correspond to GUI input.

calObjMag_48k.userEmail             = [choice.initials '@electro.dtu.dk'];
calObjMag_48k.boothName             = choice.booth;
calObjMag_48k.calculateFilterClass  = @HeaCalcFilterMag;
calObjMag_48k.channelRec            = choice.InChan;
calObjMag_48k.fs                    = choice.fs; 
calObjMag_48k.listRefSPLCalibrator  = repmat(choice.level,1,2);
calObjMag_48k.listChannelsName      = {[choice.phone ' Left'],[choice.phone ' Right']};

switch choice.phone
    case 'HD650'
    calObjMag_48k.listChannelsPlay  = [1 2];
    case 'HDA200'
    calObjMag_48k.listChannelsPlay  = [3 4];
    case 'ER2'
    calObjMag_48k.listChannelsPlay  = [5 6];
    case 'ER3'
    calObjMag_48k.listChannelsPlay  = [5 6]; 
end

%% Calibrate
calObjMag_48k.calibrate();


%% Save
save(filename,'calObjMag_48k')
save(filename_current,'calObjMag_48k')
