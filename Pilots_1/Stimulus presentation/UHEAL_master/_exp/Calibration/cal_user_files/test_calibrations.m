function test_calibrations()
% opens up a dialog for testing most recent calibrations

    path            = pwd;
    addpath(genpath([path(1:end-14)]));

    choiceTT.filename  	= 'ER3_MAG_48k_current.mat';
    choiceTT.ear        = 'l';
    choiceTT.targetlevel= 85;
    choiceTT.devID      = find_ASIO_devID;
    choiceTT.InChan     = 4;
    choiceTT.signal     = [];
    
    
    
    currentFullPath     = mfilename('fullpath');
    currentPath         = currentFullPath(1:end-length(mfilename)-1)
    if ~strcmp(currentPath,pwd),
        cd(currentPath)
    end

    calib_filenames     = dir('*.mat');
    filenames           = {calib_filenames.name};

    d                   = dialog('Position',[400 300 250 450],'Name','Test Calibration');
    
    title_txt           = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 420 210 20],...
           'String','** HEA CALTEST **','fontsize',12);
       
    date_txt            = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 400 210 20],...
           'String',date);
       
    calib_file_text     = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 340 210 40],...
           'String','Choose Calibration File');
       
    calib_file_popup    = uicontrol('Parent',d,...
           'Style','popup',...
           'Position',[45 320 160 40],...
           'String',filenames,...
           'Callback',@calib_file_callback);  
       
    ear_text         = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 280 210 40],...
           'String','Select Channel');
    ear_toggle    	= uicontrol('Parent',d,...
           'Style','togglebutton','value',0,...
           'Position',[105 270 40 30],...
           'String','L','Callback',@ear_toggle_callback);
       
    audiom_text         = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 190 210 40],...
           'String','Select Signal');
       
    audiom_push         = uicontrol('Parent',d,...
           'Style','push',...
           'Position',[75 185 100 25],...
           'String','Audiometry Freq.',...
           'Callback',@audiom_callback);   

    sweep_push          = uicontrol('Parent',d,...
            'Style','push',...
            'Position',[75 155 100 25],...
            'String','Sweep',...
            'Callback',@sweep_callback);  

     
     
    devID_txt           = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 50 210 20],...
           'String',['ASIO device ID: ' num2str(choiceTT.devID) ]);
        
   close_GUI            = uicontrol('Parent',d,...
           'Position',[89 10 70 25],...
           'String','Close',...
           'Callback',@close_GUI_callback);


    % Wait for d to close before running to completion
    uiwait(d);
   


        function calib_file_callback(calib_file_popup,callbackdata)
            idx = calib_file_popup.Value;
            choiceTT.filename  = filenames{idx};         
        end
    
        function ear_toggle_callback(ear_toggle_togglebutton,callbackdata)
            if ear_toggle_togglebutton.Value == 0,
                set(ear_toggle_togglebutton,'String','L')
                choiceTT.ear    = 'l';
            elseif ear_toggle_togglebutton.Value == 1,
                set(ear_toggle_togglebutton,'String','R')
                choiceTT.ear    = 'r';
            end
        end

        function audiom_callback(audiom_push,callbackdata)

            fs              = 48000;
            if ~isempty(strfind(choiceTT.filename, '44100'))
                fs = 44100;
            end
            frequencies     = [125 250 500 1000 2000 4000 8000];
            duration        = 5;
            overlap         = round(fs/8);
            v_t             = 0:1/fs:duration-1/fs;
            hann_temp       = hann(2*overlap);
            hanwin          = ones(fs*duration,1);

            hanwin(1:length(hann_temp)/2)          = hann_temp(1:length(hann_temp)/2);
            hanwin(1+end-length(hann_temp)/2:end)  = hann_temp(1+length(hann_temp)/2:end);
           
            for i_freq      = 1:length(frequencies)
                clc
                disp(['Press any key to play: ', num2str(frequencies(i_freq))]); 
                pause(2)
                current_signal      = hanwin.*sin(frequencies(i_freq).*2*pi*v_t');
                choiceTT.signal                 = current_signal;
                [output_signal,output_channel]	= calfilter(choiceTT.ear, choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
                % create playback object and play the sound
                player                  = HeaAudioPsych('fs', fs, 'channelsPlay', [output_channel], 'deviceID', choiceTT.devID);
                player.play(output_signal)
            end

            clc
            disp('Done!')
        end
    
        function sweep_callback(sweep_push,callbackdata)

            fs              = 48000;
            if ~isempty(strfind(choiceTT.filename, '44100'))
                fs = 44100;
            end            
            duration        = 8;
     
            overlap         = round(fs/8);
            hann_temp       = hann(2*overlap);

    
            v_t             = 0:1/fs:duration-1/fs;
            fo              = 50;
            f1              = 17000;
            sweep           = chirp(v_t',fo,max(v_t),f1,'logarithmic');
            sweep_plusMark  = [sweep; sin(1000*2*pi*[0:1/fs:3]')];
            
            hanwin          = ones(length(sweep_plusMark),1);
            hanwin(1:length(hann_temp)/2)          = hann_temp(1:length(hann_temp)/2);
            hanwin(1+end-length(hann_temp)/2:end)  = hann_temp(1+length(hann_temp)/2:end);
            
            signal          = hanwin.*sweep_plusMark;
            
           	choiceTT.signal                 = signal;
            [output_signal,output_channel]	= calfilter(choiceTT.ear, choiceTT.filename, choiceTT.targetlevel, choiceTT.signal, fs);
            % create playback object and play the sound
            player                  = HeaAudioPsych('fs', fs, 'channelsPlay', [output_channel], 'deviceID', choiceTT.devID);
            player.play(output_signal)
        end

        function close_GUI_callback(close_GUI,callbackdata)
            clc
            delete(gcf)
        end
     	
    
end
