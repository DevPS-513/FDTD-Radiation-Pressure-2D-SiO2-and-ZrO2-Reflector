function [ output_args ] = change_handle_visiblity( h_call,event_data,handle_to_change )
% This function will take the input from a specified text box
% probably on the same figure as the button that called it
% and it will load that figure assuming tab-delimited x and y values
% it will then replace the current 'XDATA' and 'YDATA' with
% the data from the file

% h is the global handle, it should have everything

% get the string

if (strcmp(get(handle_to_change,'Visible'),'on'))    
set(handle_to_change,'Visible','off')
else
set(handle_to_change,'Visible','on')
end




end

