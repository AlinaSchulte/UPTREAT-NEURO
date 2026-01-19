function measuredRMSMic = measureCalibrator(audioFunc,fs,refSPL,channelRec,...
    durationCalibrator)
%MEASURECALIBRATOR measureCalibrator records the calibrator
%using the microphone specified in 'channelRec' and returns the
%measured RMS value
%
%Function for the HeadPhone coupler

reMeasure = true;

while reMeasure
    % Make sure the setup is ready
    hWarn = warndlg('Warning: Please make sure that the sound card attenuation is set to 0 dB!','!WARNING!');
    waitfor(hWarn);
    hWarn = warndlg(sprintf('Attach calibrator to coupler and check that meter reading is %3.2f dB SPL.',refSPL),'PLACE CALIBRATOR');
    waitfor(hWarn);
    hWarn = warndlg('Make sure calibrator will stay on for at least 10 seconds before pressing OK.','PLACE CALIBRATOR');
    waitfor(hWarn);
    
    % Create the audio player 
    player = audioFunc('channelsPlay',[],'channelsRec',channelRec,'fs',fs);
    
    % Record the calibration signal
    calSig = player.rec(fs*durationCalibrator);
    
    %Cut the beginning and the end (half a second)
    calSig = calSig(round(fs/2)+1:end-round(fs/2));
    
    % Calculate root mean square value
    measuredRMSMic = rms(calSig);
    
    % Plot the signal
    figure
    plot(calSig)
    xlabel('time (samples)')
    ylabel('recorded amplitude')
    title('Calibrator signal visual check')
    
    % Re-measure calibration tone?
    button = questdlg({'Please check the calibration signal.' 'Do you want to re-calibrate the reference RMS?'},...
        'Re-record calibration tone?','Yes','No','No');
    reMeasure = strcmpi(button,'yes');
end

end