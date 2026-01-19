function [ RMSMics, dates_num ] = track_calibrations( booth, phone, from_date )
%TRACK_CALIBRATIONS reads all calibration objects from the tracker
%directory and plots measured RMS at the Mic over date. 
%
%
%   Call like: [RMSMics] = track_calibrations( 'PSY1', 'HDA200', 170404 )
%                          track_calibrations( 'PSY2', 'HD650')
%
%   Input arguments:    booth (string)
%                       phone (string)
%                       fromdate (numeric)
%
%   Output arguments:   RMSMics (NxM Matrix) 
%                               /N: number of calibrations 
%                               /M: 1(left), 2(right)
%                       dates_num (numeric, the respective dates)
%
%   Last UPDATED: 12.06.17 - Lab Manager. 

    cd('C:\data\Calibration\cal_user_files\')

    path            = pwd;
    addpath(genpath([path(1:end-14)]));

    abort_flag      = 0;
    if nargin<3,
        from_date   = 160101;
        if nargin<2,
            abort_flag  = 1;
            if nargin<1,
                warndlg('You need to select booth and phone. See the MATLAB help /doc')
            else
                warndlg('No Phone Selected. See help /doc')
            end
        end
    end
    
    if ~abort_flag

        selected_booth  = booth;
        selected_phone  = phone;
        selected_date   = from_date;

        % get filenames in calibration tracker directory.
        filenames       = dir([path(1:end-14) 'calibration_tracker']);

        RMSMics         = NaN(1,2);
        dates_num       = NaN(1);
        dates_str       = cell(1);

        ii              = 0;
        for i_f         = 1:length(filenames),

            % check if string is refering to an actual *.mat struct.
            if length(filenames(i_f).name)>3,
               current_name     = filenames(i_f).name;
               spacer_position  = findstr(current_name,'_');
               current_booth    = current_name(1:spacer_position(1)-1);

               % filter for selected booths.
               if strcmp(current_booth,selected_booth),

                    current_phone    = current_name(spacer_position(1)+1:spacer_position(2)-1);
                    %filter for selected phones.
                    if strcmp(current_phone,selected_phone),

                        current_date_str    = current_name(spacer_position(4)+1:spacer_position(5)-1);
                        current_date_num    = str2double(current_date_str);
                        % filter for lower date range
                        if current_date_num > selected_date,
                            ii      = ii+1;

                            load([path(1:end-14) 'calibration_tracker' filesep current_name])
                            
                            % Get RMSMic levels
                            
                            [RMSMics(ii,1), RMSMics(ii,2)]  = calObjMag_48k.objImpResp.measuredSPLChannel;
                            
                            % if system was NOT calibrated w 1 kHz Piston, add impulse response offset to RMSMic:
                                OneKHz          = 1000;
                                if calObjMag_48k.objImpResp(1).refFreq ~= OneKHz;
                                    [LL, RR]    = calObjMag_48k.objInvFilter.invertedImpulseResponse;
                                    [h_LL, w_L] = freqz(LL);
                                    [h_RR, w_R] = freqz(RR);
                                    
                                    offset_L    = h_LL(round((2*OneKHz)/calObjMag_48k.fs*length(w_L)));
                                    offset_R    = h_RR(round((2*OneKHz)/calObjMag_48k.fs*length(w_R)));
                                    % apply IR offset to overall RMSMic:
                                    RMSMics(ii,:)                   = RMSMics(ii,:) - 20*log10(abs([offset_L offset_R]));
                                end
                                
                            dates_str{ii}                 	= current_date_str;
                            dates_num(ii)                   = current_date_num;
                        end
                    end
                end
           end

        end

        if isnan(RMSMics(1)),
            warndlg(['No data for ' selected_booth ': ' selected_phone ', found'])
        else
            figure()
            set(gcf,'Position',[300 450 800 350])
            plot(1:length(dates_str),RMSMics,'-o','linewidth',2)
            title(['Calibration Tracking for: ' selected_booth ': ' selected_phone ])
            set(gca,'Xtick',1:length(dates_str),'XtickLabel',dates_str,'Fontsize',13)
            grid on
            ylabel('Level [db SPL] for digital RMS = 1')
            xlabel('Date calibrated')
            legend('Left','Right')
        end
    end
end
            
    
    
   
