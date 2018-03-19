function [ output_args ] = clear_row(h_event,obj_data,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


 currentrow=get(h,'Userdata');

 
 if currentrow>1
 
gui_output_table_data=get(h,'data');

[num_rows,num_cols]=size(gui_output_table_data);



gui_output_table_data(currentrow-1,:)= 0;


set(h,'data',gui_output_table_data)
set(h,'UserData',currentrow-1);
 end

end

