function [ output_args ] = write_xy_dat_from_handle(h_call,event_object,h_plot,filename,norm_on)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here




x_data=get(h_plot,'XData' );
y_data=get(h_plot,'YData' );


if norm_on==1
    y_data=y_data./(max(max(abs(y_data))));
end

write_dat(filename,x_data,y_data);
disp(['written  to: ' filename]);



end

