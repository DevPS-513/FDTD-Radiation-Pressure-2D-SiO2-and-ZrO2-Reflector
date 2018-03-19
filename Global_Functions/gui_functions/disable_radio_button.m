function [ output_args ] = disable_radio_button( h_call,event_data,h_enable,h_disable )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

val=get(h_enable,'Value');


if val==1
   
    set(h_disable,'Value',0);
    
end

end

