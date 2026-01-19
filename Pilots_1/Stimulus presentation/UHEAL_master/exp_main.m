
%% Main experimental script
clear
try
    exp_startup
    f = exp_gui;     
    uiwait(f)    
        
    exp_dorun(dat) 
           
catch
    psychrethrow(psychlasterror);
end
  