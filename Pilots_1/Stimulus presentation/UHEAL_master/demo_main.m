           
%% Main experimental script
clear
try
    exp_startup
    f = demo_gui;        
    uiwait(f)       
       
    demo_dorun(dat)
       
catch
    psychrethrow(psychlasterror);
end 
