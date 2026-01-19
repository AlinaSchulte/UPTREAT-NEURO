function [ ID ] = find_ASIO_devID( psychID )
%FIND_ASIO_DEVID checks the current audio device listings for the ASIO
%device and returns the ID. It does not require ANY input.
    InitializePsychSound
    devices     = PsychPortAudio('GetDevices');
    dev_idx     = [devices.DeviceIndex];
    dev.name    = {devices.DeviceName};
    logical_idx = cellfun(@isempty,[strfind(dev.name,'ASIO')]);
    ID          = dev_idx(~logical_idx);
    
    
    %% if find_ASIO_devID is called with a current ID, and if this ID is different to that one found, a question dialog appears.
    if nargin == 1,
        dialog      = ['Current Device is set to ' dev.name(psychID) '. The ASIO Device ID is ' num2str(ID) '. Do you want me to change it?'];
        if ID ~= psychID,
            choice = questdlg(dialog,'Device ID mismatch!','Keep','Change','Change');
            switch choice
                case 'Keep'
                    ID  = psychID;
            end
        end
    end
    
   
    
end

