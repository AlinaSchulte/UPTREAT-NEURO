function [f]=exp_gui

dat.subject=NaN;
dat.experiment_time=clock;
dat.measures=[1 2 3];
dat.measurenms = {'Pure-tone FFR' 'Click ABR' 'AEP'};
dat.ear=[2];
dat.earnms={'Left','Right'};



figpos = [0   0   600   600];
f = figure('Position',figpos,...
    'MenuBar','none','NumberTitle','off','Visible','off');

bgcol = get(gcf,'color');

uicontrol('Style','edit','units','normalized','pos',[.1 .8 .25 .1],'String','Enter ID','fontsize',15,'Callback',@idText_Callback);

uicontrol('Style','text','units','normalized','pos',[.1 .63 .25 .1],'String','Measurements:','fontsize',10,'BackgroundColor',bgcol,'HorizontalAlignment','left');
uicontrol('Style','text','units','normalized','pos',[.4 .63 .25 .1],'String','Ear:','fontsize',10,'BackgroundColor',bgcol,'HorizontalAlignment','left');
uicontrol('Style','listbox','String',dat.measurenms,'Max',100,'Min',1,...
    'Value',[1 2 3],'units','normalized','pos',[.1 .6 .25 .1],'Callback',@condition_Callback);
uicontrol('Style','listbox','String',dat.earnms,'Max',3,'Min',1,...
    'Value',[2],'units','normalized','pos',[.4 .6 .25 .1],'Callback',@earlist_Callback);

uicontrol('Style','pushbutton','units','normalized','pos',[.4 .8 .25 .1],'fontsize',14,'String','RUN TEST','Callback',@runbutton_Callback);

movegui(f,'center')
set(f,'Visible','on')

%%
    function condition_Callback(hObject,EventData)
        id_selected = get(hObject,'Value');
        %list = get(hObject,'String');
        %dat.fcs=str2double(list(1:end-1))';
        dat.measures=id_selected;
    end



    function earlist_Callback(hObject,EventData)
        dat.ear = get(hObject,'Value');
    end

    function idText_Callback(hObject,EventData)
        dat.subject = get(hObject,'string');
    end


    function runbutton_Callback(hObject,EventData)
        if isnan(dat.subject)
            tmp = inputdlg('Subject ID:');
            dat.subject=tmp{1};
        end
        assignin('base','dat',dat)
        close(f)
    end

end