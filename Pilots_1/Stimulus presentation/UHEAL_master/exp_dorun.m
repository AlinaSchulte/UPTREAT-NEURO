function exp_dorun(dat)

if ismember(1,dat.measures)
ffr_SW_main(dat) %trigger values 10-40
%pt_ffr_main(dat)
pause
end
if ismember(2,dat.measures)
click_abr_main_2xramp(dat) %trigger values 90-120
pause
end
if ismember(3,dat.measures)
AEP_EFR_main_4ISI(dat) %trigger values 130
pause
end

end