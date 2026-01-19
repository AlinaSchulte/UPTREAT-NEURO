%% functions
function saveFFRdata(data,savedir,subid)


    %%  Save mat
    savefile = [savedir filesep subid '_FFR_4Hz.mat'];
    save(savefile,'data','-v7.3');

end
