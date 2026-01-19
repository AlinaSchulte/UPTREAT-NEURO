%% script extracting iphone signals and the template sgnals (before playing on iphone) and calculating IR

cd 'C:\data\Calibration'
addpath(genpath(pwd))

load template_refsig
load template_sweep
N = length(template_sweep);

% get templates

for idR = 1:3
    
    % sweep sig             
    eval(['load(''sweep_iPhone_L', num2str(idR), '.mat'')'])
    [xc, lag] = xcorr(signal, template_sweep);
    [xcv, ids_max] = sort(xc, 'descend');
    vmax_5(1) = xcv(1);    
    idsmax_5(1) = ids_max(1);
    lags_5(1) = lag(idsmax_5(1));
    sweeps_5_3(1:N, 1, idR) = signal(lags_5(1) + [1:length(template_sweep)]);
    for id = 2:5
        rem_ids = ids_max > idsmax_5(id-1) - 10000 & ids_max < idsmax_5(id-1) + 10000;
        ids_max = ids_max(~rem_ids);
        xcv = xcv(~rem_ids);
        vmax_5(id) = xcv(1);
        idsmax_5(id) = find(xc == vmax_5(id));
        lags_5(id) = lag(idsmax_5(id));
        sweeps_5_3(1:N, id, idR) = signal(lags_5(id) + [1:length(template_sweep)]);
    end

    % refSig
    eval(['load(''refSig_iPhone_L', num2str(idR), '.mat'')'])
    [xc, lag] = xcorr(signal, template_refsig);
    [xcv, ids_max] = sort(xc, 'descend');
    vmax(1) = xcv(1);    
    idsmax(1) = ids_max(1);
    lags(1) = lag(idsmax(1));
    refsig_3(1:length(template_refsig), idR) = signal(lags(1) + [1:length(template_refsig)]);    
end

%%

load obj

% Create reference signal
refSig = obj.refAmp*sin((0:1:round(obj.durationRecording*obj.fs)-1)/obj.fs*2*pi*obj.refFreq)';
% Measure headphone frequency response
HPSig = zeros((obj.durationRecording-1)*obj.fs,obj.nPlace);
rmsHPSig = zeros(1,obj.nPlace);

H = zeros(ceil((obj.preDelay + obj.decay)/1000*obj.fs/2),obj.nPlace);
h = zeros(ceil((obj.preDelay + obj.decay)/1000*obj.fs),obj.nPlace);


for ll = 1:3
    
    clc; fprintf('Measuring Channel %s \n',obj.channelPlay);
    
    % Play and record reference signal for "obj.nPlace" times
    Sig = refsig_3(:, ll);    
    HPSig(:,ll) = Sig(obj.fs/2+1:(obj.durationRecording-1)*obj.fs+obj.fs/2);
    
    % Calculate root mean square value of recorded signal
    rmsHPSig(:,ll) = rms(HPSig(:,ll));
    
    % Measure impulse response
    [IR(:,ll),x,y,t,obj.nRepSweep,D,obj.fs] = measIR_mfe(sweeps_5_3(:,:,idR), obj, obj.nRepSweep,obj.sweepLen,obj.fs,obj.scale);
    
    %Produces a one-sided amplitude spectrum from an impulse response
    [H(:,ll),h(:,ll),fNorm] = getH_mfe(obj, obj.preDelay,obj.decay,IR(:,ll),obj.fs,obj.refFreq);    
   
end

rmsHPSigmean = mean(rmsHPSig);
obj.measuredSPLChannel = obj.refSPLCalibrator+20*log10(sqrt(2)*rmsHPSigmean/(obj.refAmp*obj.measuredRMSMic));
disp(obj.measuredSPLChannel)

% Calculate mean values of TFs for calibrating Channel
HCh = mean(abs(H(:,1:obj.nPlace)),2);      % mean (measured) Magnitude of left EP
obj.measuredImpulseResponse = mean(h(:,1:obj.nPlace),2);

% save the object here eg as obj_post and substitute it in HeaCalPSYBooth  just before
% executing line 150 where inverse filter is created

                