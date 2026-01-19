function [ channelID, ID] = find_INPUT_channel()
%FIND_INPUT_CHANNEL returns ID (device) and channel ID (Input).
%   Call like:  [ channelID, devID ] = find_INPUT_channel()
%
%   s156019@student.dtu.dk
%
    addpath(genpath(pwd))
    InitializePsychSound

    % GET ASIO DEVICE ID
    devices     = PsychPortAudio('GetDevices');
    dev_idx     = [devices.DeviceIndex];
    dev.name    = {devices.DeviceName};
    logical_idx = cellfun(@isempty,[strfind(dev.name,'ASIO')]);
    ID          = dev_idx(~logical_idx);
    ID = find_ASIO_devID;

    % GET NUMBER OF DEVICE INPUT CHANNELS
    channel_num = devices(~logical_idx).NrInputChannels;

    % EIGEN-PREFERENCES
    maxsecs     = 5;   % duration of measurement
    fs          = 48000;


    search      = true;

        while search
        close all
        % INITIALIZE
        hWarn = warndlg('Make some noise at the microphone after pressing OK','PLACE CALIBRATOR / MAKE NOISE');
        waitfor(hWarn);
        pahandle = PsychPortAudio('Open', ID, 2, 0, fs, channel_num);
        PsychPortAudio('GetAudioData', pahandle, 10);
        PsychPortAudio('Start', pahandle, 0, 0, 1);
        fprintf('Audio capture started. Duration: %i sec.\n',maxsecs);
        % We retrieve status once to get access to SampleRate:
        s = PsychPortAudio('GetStatus', pahandle);
        recordedaudio = [];

        while ((length(recordedaudio) / s.SampleRate) < maxsecs)
            % Wait a second...
            WaitSecs(1);

            % Query current capture status and print it to the Matlab window:
            s = PsychPortAudio('GetStatus', pahandle);

            % Print it:
            fprintf('Audio capture started. Time left: %i sec.\n',round(maxsecs-(length(recordedaudio) / s.SampleRate)));

            % Retrieve pending audio data from the drivers internal ringbuffer:
            audiodata = PsychPortAudio('GetAudioData', pahandle);
            nrsamples = size(audiodata, 2);

            % And attach it to full sound vector:
            recordedaudio = [recordedaudio audiodata]; %#ok<AGROW>
        end

        % Stop capture:
        PsychPortAudio('Stop', pahandle);

        % Perform a last fetch operation to get all remaining data from the capture engine:
        audiodata = PsychPortAudio('GetAudioData', pahandle);

        % Attach it to our full sound vector:
        recordedaudio = [recordedaudio audiodata];

        % Close the audio device:
        PsychPortAudio('Close', pahandle);

        audiosum        = sum(abs(audiodata),2);
        channelID       = find(audiosum==max(audiosum));

        checkfig    = figure();
        plot(audiodata(channelID,:))
        xlabel('time/samples')
        ylabel('Amplitude')

        qststr  = ['Input signal for Device ID: ' num2str(ID) ' at channel: ' num2str(channelID)];
        title(qststr)
        choice = questdlg(qststr,'Check Input','Repeat Search','Start Calibration','Start Calibration');
        % Handle response
        
        switch choice
            case 'Start Calibration'
                close(checkfig)
                answ    = 'n';
                
            case 'Repeat Search'
                disp([choice ' coming right up.'])
                answ = 'y';
        end

            search = any(strcmpi(answ,{'yes' 'y'}));
    end

end

