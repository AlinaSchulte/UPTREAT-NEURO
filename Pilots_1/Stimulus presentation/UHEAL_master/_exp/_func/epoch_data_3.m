function [data_out] = epoch_data_3(data_in,trial_id,nr_epochs)
%Epoch data

%More info following

if ~isstruct(data_in)
    error('Data not in right format');
else
    clc
    disp('Epoching data');
end
if length(trial_id)>1
    idx = find(data_in.trialinfo==trial_id(1) | data_in.trialinfo==trial_id(2));
else
    idx=find(data_in.trialinfo==trial_id);
end
if length(nr_epochs)==1
    idx_tmp=idx(1:nr_epochs);
else
    idx_tmp=nr_epochs;
end

% generate epochs
epoched_data = [];
chans = size(data_in.trial{1},1);
times = size(data_in.trial{1},2);
epoched_data = reshape(cell2mat(data_in.trial(idx_tmp)),chans,times,length(idx_tmp)); % epoch x chan x time
clc
disp(['trials done.'])


%if size(epoched_data,1)>1 %mean over channels
%    epoched_data = mean(epoched_data,2);
%end
data_out = permute(epoched_data,[3,1,2]);

end
