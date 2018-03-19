function [ output_args ] = plot_from_visible_textbox( h_call,event_data,h_textbox,h_plot )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if(strcmp(get(h_textbox,'Visible'),'on'))

path=get(h_textbox,'string');

if exist(path)
data=load(path);
set(h_plot,'XDATA',data(:,1));
set(h_plot,'YDATA',data(:,2));

set(h_textbox,'Visible','off');
end

end


set(h_textbox,'string','{no ./path}');
end

