classdef HEATriggerbox < handle & matlab.mixin.CustomDisplay
    
    properties (Access = private)
        com_port_name
        com_port_handle
        com_port_isopen
    end
    
    methods
        
        function obj = HEATriggerbox()
            obj.com_port_name = '';
            obj.com_port_handle = 0;
            obj.com_port_isopen = 0;
        end
        
        
        function find_triggerbox_win(obj)

            disp('Searching for devices ...');

            [~,res] = system('wmic path Win32_PnPEntity get DeviceID, Caption /format:csv');
            rows = textscan(res, '%s', 'Delimiter', '\n');
            rows = rows{1};
            for row_i = 1:length(rows)
                row = rows{row_i};
                if isempty(row)
                    continue
                end

                scan = textscan(row, '%s', 'Delimiter', ',');
                fields = scan{1};
                start_i = row_i;
                break;
            end

            devid_field_idx = find(cellfun(@(x) strcmp(x, 'DeviceID'), fields));
            caption_field_idx = find(cellfun(@(x) strcmp(x, 'Caption'), fields));

            found = 0;
            for row_i = start_i:length(rows)
                row = rows{row_i};
                scan = textscan(row, '%s', 'Delimiter',',');
                values = scan{1};

                if isempty(strfind(values{devid_field_idx}, 'VID_16C0')) && ...
                        isempty(strfind(values{devid_field_idx}, 'PID_0483'))
                    continue
                end

                fprintf('Found a device with a suitable ID: %s\n', values{caption_field_idx});
                caption = values{caption_field_idx};
                found = 1;
                break; 
            end

            if found == 0
                disp('Could not find the triggerbox. Check the connection');
            else
                disp('---> looking for a matching serial device ...');

                [~,res] = system('wmic path Win32_SerialPort get DeviceID, Caption /format:csv');
                rows = textscan(res, '%s', 'Delimiter', '\n');
                rows = rows{1};
                for row_i = 1:length(rows)
                    row = rows{row_i};
                    if isempty(row)
                        continue
                    end

                    scan = textscan(row, '%s', 'Delimiter', ',');
                    fields = scan{1};
                    start_i = row_i;
                    break;
                end
                devid_field_idx = find(cellfun(@(x) strcmp(x, 'DeviceID'), fields));
                caption_field_idx = find(cellfun(@(x) strcmp(x, 'Caption'), fields));

                for row_i = start_i:length(rows)
                    row = rows{row_i};
                    scan = textscan(row, '%s', 'Delimiter',',');
                    values = scan{1};

                    if strcmp(caption, values{caption_field_idx})
                        fprintf('Match found. Setting port to %s\n', values{devid_field_idx});
                        obj.com_port_name = values{devid_field_idx};
                        break;
                    end
                end
            end
            
        end
        
        
        function success = connect(obj, port)
            
            if nargin == 2
                obj.com_port_name = port;
            elseif nargin == 1
                if isempty(obj.com_port_name)
                    fprintf('No port name given. Please provide a COM port name or run the method "find_triggerbox_win()"\n');
                    return
                end
            end
            
            try
                obj.com_port_handle = serial(obj.com_port_name, 'BaudRate', 9600);
                fopen(obj.com_port_handle);
                obj.com_port_isopen = 1;
                fprintf('Triggerbox at %s connected\n', obj.com_port_name);
                success = 1;
            catch
                obj.com_port_is_open = 0;
                fprintf('Could not connect to triggerbox at %s\n', obj.com_port_name);
                success = 0;
            end
            
        end
        
        
        function connected = is_connected(obj)
            connected = obj.com_port_isopen;
        end
        
        
        function latencies = estimate_latency(obj, repetitions, interval)
            if nargin < 2
                repetitions = 20;
                interval = 0.5;
            end
            triggercode = 12;
            
            fprintf('Updating trigger %d times, this will take about %.1f seconds\n', repetitions, repetitions*interval);
            latencies = zeros(repetitions, 1);
            if obj.com_port_isopen
                for j=1:repetitions
                    tstart = tic;
                    out = obj.set_trigger(triggercode);
                    assert(strcmp(out, num2str(triggercode)));
                    telapsed = toc(tstart);
                    pause(interval);
                    latencies(j) = telapsed;
                end
                fprintf('\n');
                
                fprintf('It takes at worst %.3f ms (avg: %.3f ms) to update a trigger\n', 1000*max(latencies), 1000*mean(latencies));
            else
                fprintf('The triggerbox has not been connected. Run the method "connect_triggerbox([port])" first\n');
            end
        end
        
        
        function out = set_trigger(obj, trigger_value)
            if obj.com_port_isopen
                fwrite(obj.com_port_handle, trigger_value, 'uint8');
                while obj.com_port_handle.BytesAvailable == 0
                end
                out = char(fread(obj.com_port_handle, obj.com_port_handle.BytesAvailable)');
            else
                ME = MException('HEATriggerbox:notConnected', ...
                    'No connection to triggerbox, cannot set trigger.');
                throw(ME);
            end
            
        end
        
        
        function disconnect(obj)
            if obj.com_port_isopen
                fclose(obj.com_port_handle);
                obj.com_port_isopen = 0;
            end
        end
        
        
        function delete(obj)
            if obj.com_port_isopen
                fclose(obj.com_port_handle);
                obj.com_port_isopen = 0;
            end
        end
        
        
        function status(obj)
            disp(obj);
        end
        
    end
    
    
    methods (Access = protected)
        
        function header = getHeader(obj)
            if ~isscalar(obj)
                header = getHeader@matlab.mixin.CustomDisplay(obj);
            else
                headerStr = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
                if isempty(obj.com_port_name)
                    headerStr = [headerStr, ' with no serial port set.'];
                else
                    headerStr = [headerStr, ' at ', obj.com_port_name, '.'];
                end
                
                if obj.com_port_isopen
                    headerStr = [headerStr, ' Connected.'];
                else
                    headerStr = [headerStr, ' Not connected.'];
                end
                
                header = sprintf('%s\n', headerStr);
                
            end

        end
        
    end
    
end