 clearvars all; close all; clc

global choice

addpath(genpath(pwd))
cd('C:\data\jonmarc\UHEAL_master\_exp\Calibration')
startupHeaAudio
startupCalibration

choice          = dialogRR();
choice.path     = pwd;
choice.fs = 48000;
if isfield(choice,'calib')
    calScript

    close all;
    disp(['Calibrated Phone: ' choice.phone '. Data saved as current file and in tracker folder.'])
    test_calibrations
else
    clc
end
