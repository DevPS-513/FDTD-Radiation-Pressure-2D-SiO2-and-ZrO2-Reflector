function [ output_args ] = copy_handle( h_call,event_datam,handle_to_copy,handle,parent)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

new_handle=copyobj(handle_to_copy,parent);

set(new_handle,'Visible','on')

end

